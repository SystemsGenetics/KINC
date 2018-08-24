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



def load_cmx(filename, num_genes, num_clusters):
	netlist = pd.read_csv(args.INPUT, sep="\t", header=None)
	cmx = np.zeros((num_genes * num_clusters, num_genes * num_clusters), dtype=np.float32)

	for idx in xrange(len(netlist.index)):
		i = netlist.iloc[idx, 0]
		j = netlist.iloc[idx, 1]
		k = netlist.iloc[idx, 2]
		r = netlist.iloc[idx, 9]

		cmx[i * num_clusters + k, j * num_clusters + k] = r
		cmx[j * num_clusters + k, i * num_clusters + k] = r

	return cmx



def compute_pruned_matrix(S, threshold):
	S_pruned = np.copy(S)
	S_pruned[abs(S) < threshold] = 0
	S_pruned = S_pruned[~np.all(S_pruned == 0, axis=1)]
	S_pruned = S_pruned[:, ~np.all(S_pruned == 0, axis=0)]

	return S_pruned



def powerlaw(args):
	# load correlation matrix
	S = load_cmx(args.INPUT, args.NUM_GENES, args.MAX_CLUSTERS)

	# iterate until network is sufficiently scale-free
	threshold = args.TSTART

	while True:
		# compute thresholded adjacency matrix
		A = (abs(S) >= threshold)

		# compute degree of each node
		for i in xrange(A.shape[0]):
			A[i, i] = 0

		degrees = np.array([sum(A[i]) for i in xrange(A.shape[0])])

		# compute degree distribution
		bins = max(5, degrees.max())
		hist, _ = np.histogram(degrees, bins=bins, range=(1, bins))
		bin_edges = range(1, len(hist) + 1)

		# modify histogram values to work with loglog plot
		hist += 1

		# plot degree distribution
		if args.VISUALIZE:
			plt.subplots(1, 2, figsize=(10, 5))
			plt.subplot(121)
			plt.plot(bin_edges, hist, "ko")
			plt.subplot(122)
			plt.loglog(bin_edges, hist, "ko")
			plt.savefig("plots/powerlaw/%03d.png" % (int(threshold * 100)))
			plt.close()

		# compute correlation
		x = np.log(bin_edges)
		y = np.log(hist)

		r, p = scipy.stats.pearsonr(x, y)

		# output results of threshold test
		print("%g\t%g\t%g" % (threshold, r, p))

		# break if power law is satisfied
		if r < 0 and p < 1e-20:
			break

		# decrement threshold and fail if minimum threshold is reached
		threshold -= args.TSTEP
		if threshold < args.TSTOP:
			print("error: could not find an adequate threshold above stopping threshold")
			sys.exit(0)

	return threshold



def compute_degenerate(eigens):
	unique = []
	
	for i in xrange(len(eigens)):
		if len(unique) == 0 or abs(eigens[i] - unique[-1]) > 1e-6:
			unique.append(eigens[i])

	return unique



def compute_spacings(eigens, pace):
	# extract eigenvalues for spline based on pace
	x = eigens[::pace]
	y = np.linspace(0, 1, len(x))

	# compute spline
	spl = scipy.interpolate.splrep(x, y)

	# extract interpolated eigenvalues from spline
	spline_eigens = scipy.interpolate.splev(eigens, spl)

	# compute spacings between interpolated eigenvalues
	spacings = np.empty(len(eigens) - 1)
	
	for i in xrange(len(spacings)):
		spacings[i] = (spline_eigens[i + 1] - spline_eigens[i]) * len(eigens)

	return spacings



def compute_chi_square_pace(eigens, pace):
	# compute eigenvalue spacings
	spacings = compute_spacings(eigens, pace)

	# compute nearest-neighbor spacing distribution
	hist_min = 0.0
	hist_max = 3.0
	num_bins = 60
	bin_width = (hist_max - hist_min) / num_bins

	hist, _ = np.histogram(spacings, num_bins, (hist_min, hist_max))
	
	# compote chi-square value from nnsd
	chi = 0
	
	for i in xrange(len(hist)):
		# compute O_i, the number of elements in bin i
		O_i = hist[i]

		# compute E_i, the expected value of Poisson distribution for bin i
		E_i = (math.exp(-i * bin_width) - math.exp(-(i + 1) * bin_width)) * len(eigens)

		# update chi-square value based on difference between O_i and E_i
		chi += (O_i - E_i) * (O_i - E_i) / E_i

	print("pace: %d, chi: %g" % (pace, chi))

	return chi



def compute_chi_square(eigens):
	# compute unique eigenvalues
	unique = compute_degenerate(eigens)

	print("eigenvalues: %d" % len(eigens))
	print("unique eigenvalues: %d" % len(unique))

	# make sure there are enough eigenvalues
	if len(unique) < 50:
		return -1

	# perform several chi-square tests by varying the pace
	chi = 0
	num_tests = 0

	for pace in xrange(10, 41):
		# make sure there are enough eigenvalues for pace
		if len(unique) / pace < 5:
			break

		chi += compute_chi_square_pace(unique, pace)
		num_tests += 1
	
	# compute average of chi-square tests
	chi /= num_tests

	# return chi value
	return chi



def rmt(args):
	# load correlation matrix
	S = load_cmx(args.INPUT, args.NUM_GENES, args.MAX_CLUSTERS)

	# iterate until chi value goes below 99.607 then above 200
	final_threshold = 0
	final_chi = float("inf")
	max_chi = -float("inf")
	threshold = args.TSTART

	while max_chi < 200:
		# compute pruned matrix
		S_pruned = compute_pruned_matrix(S, threshold)

		# make sure pruned matrix is not empty
		chi = -1

		if S_pruned.shape[0] > 0:
			# compute eigenvalues of pruned matrix
			eigens, _ = np.linalg.eigh(S_pruned)

			# compute chi-square value from NNSD of eigenvalues
			chi = compute_chi_square(eigens)

		# make sure chi-square test succeeded
		if chi != -1:
			# save most recent chi-square value less than critical value
			if chi < 99.607:
				final_chi = chi
				final_threshold = threshold

			# save largest chi-square value which occurs after final_chi
			if final_chi < 99.607 and chi > final_chi:
				max_chi = chi

		# output results of threshold test
		print("%f\t%d\t%f" % (threshold, S_pruned.shape[0], chi))

		# decrement threshold and fail if minimum threshold is reached
		threshold -= args.TSTEP
		if threshold < args.TSTOP:
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
	parser.add_argument("-i", "--input", required=True, help="correlation matrix file", dest="INPUT")
	parser.add_argument("--genes", type=int, required=True, help="number of genes", dest="NUM_GENES")
	parser.add_argument("--method", default="rmt", choices=["powerlaw", "rmt"], help="thresholding method", dest="METHOD")
	parser.add_argument("--tstart", type=float, default=0.99, help="starting threshold", dest="TSTART")
	parser.add_argument("--tstep", type=float, default=0.001, help="threshold step size", dest="TSTEP")
	parser.add_argument("--tstop", type=float, default=0.5, help="stopping threshold", dest="TSTOP")
	parser.add_argument("--minclus", type=int, default=1, help="minimum clusters", dest="MIN_CLUSTERS")
	parser.add_argument("--maxclus", type=int, default=5, help="maximum clusters", dest="MAX_CLUSTERS")
	parser.add_argument("--visualize", action="store_true", help="whether to visualize results", dest="VISUALIZE")

	args = parser.parse_args()

	# print arguments
	pprint.pprint(vars(args))

	# load data
	cmx = pd.read_csv(args.INPUT, sep="\t")

	# initialize method
	compute_threshold = METHODS[args.METHOD]

	print(compute_threshold(args))
