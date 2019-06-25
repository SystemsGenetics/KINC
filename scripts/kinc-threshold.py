import argparse
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pprint
import scipy.interpolate
import scipy.stats
import seaborn as sns
import sklearn.cluster
import sklearn.mixture
import sys



def load_cmx(filename, n_genes, n_clusters):
	netlist = pd.read_csv(filename, sep="\t", header=None)
	cmx = np.zeros((n_genes * n_clusters, n_genes * n_clusters), dtype=np.float32)

	for idx in range(len(netlist.index)):
		i = netlist.iloc[idx, 0]
		j = netlist.iloc[idx, 1]
		k = netlist.iloc[idx, 2]
		r = netlist.iloc[idx, 9]

		cmx[i * n_clusters + k, j * n_clusters + k] = r
		cmx[j * n_clusters + k, i * n_clusters + k] = r

	return cmx



def powerlaw(args):
	# load correlation matrix
	S = load_cmx(args.input, args.n_genes, args.maxclus)

	# iterate until network is sufficiently scale-free
	threshold = args.tstart

	while True:
		# compute thresholded adjacency matrix
		A = (abs(S) >= threshold)

		# compute degree of each node
		for i in range(A.shape[0]):
			A[i, i] = 0

		degrees = np.array([sum(A[i]) for i in range(A.shape[0])])

		# compute degree distribution
		bins = max(5, degrees.max())
		hist, _ = np.histogram(degrees, bins=bins, range=(1, bins))
		bin_edges = range(1, len(hist) + 1)

		# modify histogram values to work with loglog plot
		hist += 1

		# compute correlation
		x = np.log(bin_edges)
		y = np.log(hist)

		r, p = scipy.stats.pearsonr(x, y)

		# plot degree distribution
		if args.visualize:
			plt.subplots(1, 2, figsize=(10, 5))
			plt.subplot(121)
			plt.title("Degree Distribution")
			plt.plot(bin_edges, hist, "o")
			plt.subplot(122)
			plt.title("Degree Distribution (log-log)")
			sns.regplot(x, y)
			plt.savefig("powerlaw_%03d.png" % (int(threshold * 1000)))
			plt.close()

		# output results of threshold test
		print("%0.3f\t%0.3f\t%e" % (threshold, r, p))

		# break if power law is satisfied
		if r < 0 and p < 1e-20:
			break

		# decrement threshold and fail if minimum threshold is reached
		threshold -= args.tstep
		if threshold < args.tstop:
			print("error: could not find an adequate threshold above stopping threshold")
			sys.exit(0)

	return threshold



def compute_pruned_matrix(S, threshold):
	S_pruned = np.copy(S)
	S_pruned[abs(S) < threshold] = 0
	S_pruned = S_pruned[~np.all(S_pruned == 0, axis=1)]
	S_pruned = S_pruned[:, ~np.all(S_pruned == 0, axis=0)]

	return S_pruned



def compute_unique(values):
	unique = []
	
	for i in range(len(values)):
		if len(unique) == 0 or abs(values[i] - unique[-1]) > 1e-6:
			unique.append(values[i])

	return unique



def compute_spline(values, pace):
	# extract values for spline based on pace
	x = values[::pace]
	y = np.linspace(0, 1, len(x))

	# compute spline
	spl = scipy.interpolate.splrep(x, y)

	# extract interpolated eigenvalues from spline
	spline_values = scipy.interpolate.splev(values, spl)

	return spline_values



def compute_spacings(values):
	spacings = np.empty(len(values) - 1)
	
	for i in range(len(spacings)):
		spacings[i] = (values[i + 1] - values[i]) * len(values)

	return spacings



def compute_chi_square_helper(values):
	# compute eigenvalue spacings
	spacings = compute_spacings(values)

	# compute nearest-neighbor spacing distribution
	hist_min = 0.0
	hist_max = 3.0
	num_bins = 60
	bin_width = (hist_max - hist_min) / num_bins

	hist, _ = np.histogram(spacings, num_bins, (hist_min, hist_max))
	
	# compote chi-square value from nnsd
	chi = 0
	
	for i in range(len(hist)):
		# compute O_i, the number of elements in bin i
		O_i = hist[i]

		# compute E_i, the expected value of Poisson distribution for bin i
		E_i = (math.exp(-i * bin_width) - math.exp(-(i + 1) * bin_width)) * len(values)

		# update chi-square value based on difference between O_i and E_i
		chi += (O_i - E_i) * (O_i - E_i) / E_i

	return chi



def compute_chi_square(eigens, spline=True):
	# make sure there are enough eigenvalues
	if len(eigens) < 50:
		return -1

	# use spline interpolation if specified
	if spline:
		# perform several chi-square tests with spline interpolation by varying the pace
		chi = 0
		num_tests = 0

		for pace in range(10, 41):
			# make sure there are enough eigenvalues for pace
			if len(eigens) / pace < 5:
				break

			# compute spline-interpolated eigenvalues
			eigens = compute_spline(eigens, pace)

			# compute chi-squared value
			chi_pace = compute_chi_square_helper(eigens)

			print("pace: %d, chi-squared: %g" % (pace, chi_pace))

			# compute chi-squared value
			chi += chi_pace
			num_tests += 1

		# return average of chi-square tests
		return chi / num_tests

	# perform a single chi-squared test without spline interpolation
	else:
		return compute_chi_square_helper(eigens)



def rmt(args):
	# load correlation matrix
	S = load_cmx(args.input, args.n_genes, args.maxclus)

	# iterate until chi value goes below 99.607 then above 200
	final_threshold = 0
	final_chi = float("inf")
	max_chi = -float("inf")
	threshold = args.tstart

	while max_chi < 200:
		# compute pruned matrix
		S_pruned = compute_pruned_matrix(S, threshold)

		# make sure pruned matrix is not empty
		chi = -1

		if S_pruned.shape[0] > 0:
			# compute eigenvalues of pruned matrix
			eigens, _ = np.linalg.eigh(S_pruned)

			print("eigenvalues: %d" % len(eigens))

			# compute unique eigenvalues
			eigens = compute_unique(eigens)

			print("unique eigenvalues: %d" % len(eigens))

			# compute chi-square value from NNSD of eigenvalues
			chi = compute_chi_square(eigens, spline=args.spline)

		# make sure chi-square test succeeded
		if chi != -1:
			# plot eigenvalue distribution
			if args.visualize:
				plt.subplots(1, 2, figsize=(10, 5))
				plt.subplot(121)
				plt.title("Eigenvalues")
				plt.plot(eigens, ".")
				plt.subplot(122)
				plt.title("Eigenvalue Spacing Distribution")
				plt.hist(compute_spacings(eigens), bins=60, range=(0, 3))
				plt.savefig("rmt_%03d.png" % (int(threshold * 1000)))
				plt.close()

			# save most recent chi-square value less than critical value
			if chi < 99.607:
				final_chi = chi
				final_threshold = threshold

			# save largest chi-square value which occurs after final_chi
			if final_chi < 99.607 and chi > final_chi:
				max_chi = chi

		# output results of threshold test
		print("%0.3f\t%d\t%g" % (threshold, S_pruned.shape[0], chi))

		# decrement threshold and fail if minimum threshold is reached
		threshold -= args.tstep
		if threshold < args.tstop:
			print("error: could not find an adequate threshold above stopping threshold")
			sys.exit(0)

	return final_threshold



if __name__ == "__main__":
	# define threshold methods
	METHODS = {
		"powerlaw": powerlaw,
		"rmt": rmt
	}

	# parse command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--input", help="correlation matrix file", required=True)
	parser.add_argument("--n-genes", help="number of genes", type=int, required=True)
	parser.add_argument("--method", help="thresholding method", default="rmt", choices=["powerlaw", "rmt"])
	parser.add_argument("--tstart", help="starting threshold", type=float, default=0.99)
	parser.add_argument("--tstep", help="threshold step size", type=float, default=0.001)
	parser.add_argument("--tstop", help="stopping threshold", type=float, default=0.5)
	parser.add_argument("--spline", help="whether to use spline interpolation", action="store_true")
	parser.add_argument("--minclus", help="minimum clusters", type=int, default=1)
	parser.add_argument("--maxclus", help="maximum clusters", type=int, default=5)
	parser.add_argument("--visualize", help="whether to visualize results", action="store_true")

	args = parser.parse_args()

	# print arguments
	pprint.pprint(vars(args))

	# load data
	cmx = pd.read_csv(args.input, sep="\t")

	# initialize method
	compute_threshold = METHODS[args.method]

	# compute threshold
	threshold = compute_threshold(args)

	print("%0.3f" % (threshold))
