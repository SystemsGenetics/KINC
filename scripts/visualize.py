import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy.stats
import seaborn as sns



def plot_clusdist(netlist, output_dir):
	max_clusters = netlist["Num_Clusters"].max()

	sns.distplot(netlist["Num_Clusters"], kde=False)
	plt.title("Cluster Size Distribution")
	plt.xlabel("Cluster Size")
	plt.ylabel("Frequency")
	plt.xticks(range(1, max_clusters + 1))
	plt.savefig("%s/clusdist.png" % output_dir)
	plt.close()



def plot_corrdist(netlist, output_dir):
	sns.distplot(netlist["sc"], bins=np.arange(-1, 1.01, 0.01), kde=False)
	plt.title("Correlation Distribution")
	plt.xlabel("Correlation")
	plt.ylabel("Frequency")
	plt.savefig("%s/corrdist.png" % output_dir)
	plt.close()



def plot_coverage(netlist, output_dir):
	thresholds = np.arange(0.5, 1.0, 0.01)
	coverage = np.zeros(len(thresholds))

	for i, threshold in enumerate(thresholds):
		edges = netlist[abs(netlist["sc"]) >= threshold]
		genes = set(edges["Source"]).union(set(edges["Target"]))
		coverage[i] = len(genes)

	plt.plot(thresholds, coverage)
	plt.title("Gene Coverage")
	plt.xlabel("Abs. Correlation Threshold")
	plt.ylabel("Genes")
	plt.ylim(bottom=0)
	plt.savefig("%s/coverage.png" % output_dir)
	plt.close()



def plot_pairwise(emx, netlist, output_dir, limits=None, range_args=None):
	# determine range of plots to render
	if range_args != None:
		start, stop, step = range_args
		indices = netlist.index[start:stop:step]
	else:
		indices = netlist.index

	# iterate through each network edge
	for idx in indices:
		edge = netlist.iloc[idx]
		x = edge["Source"]
		y = edge["Target"]
		k = edge["Cluster"]

		print("%-20s %-20s %d" % (x, y, k))

		# extract pairwise data
		labels = np.array([int(s) for s in edge["Samples"]])
		mask1 = (labels != 9)

		labels = labels[mask1]
		mask2 = (labels == 1)

		data = emx.loc[[x, y]].values[:, mask1]

		# highlight samples in the edge
		colors = np.array(["k" for _ in labels])
		colors[mask2] = "r"

		# compute Spearman correlation
		r, p = scipy.stats.spearmanr(data[0, mask2], data[1, mask2])

		# create figure
		plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10, 5))

		# create density plot
		plt.subplot(121)

		if limits != None:
			plt.xlim(limits)
			plt.ylim(limits)

		sns.kdeplot(data[0], data[1], shade=True, shade_lowest=False)

		# create scatter plot
		plt.subplot(122)
		plt.title("k=%d, samples=%d, spearmanr=%0.2f" % (k, edge["Cluster_Samples"], r))

		if limits != None:
			plt.xlim(limits)
			plt.ylim(limits)

		plt.xlabel(x)
		plt.ylabel(y)
		plt.scatter(data[0], data[1], color="w", edgecolors=colors)

		# save plot to file
		plt.savefig("%s/%s_%s_%d.png" % (output_dir, x, y, k))
		plt.close()



if __name__ ==  "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--emx", help="expression matrix file", required=True)
	parser.add_argument("--netlist", help="netlist file", required=True)
	parser.add_argument("--output-dir", help="output directory", default=".")
	parser.add_argument("--clusdist", help="plot cluster distribution", action="store_true")
	parser.add_argument("--corrdist", help="plot correlation distribution", action="store_true")
	parser.add_argument("--coverage", help="plot gene coverage vs correlation", action="store_true")
	parser.add_argument("--pairwise", help="plot pairwise plots", action="store_true")
	parser.add_argument("--pw-range", help="range of pairwise plots to render", nargs=3, type=int, metavar=("START", "STOP", "STEP"))
	parser.add_argument("--pw-scale", help="use the same limits for all pairwise plots", action="store_true")

	args = parser.parse_args()

	# load input data
	emx = pd.read_csv(args.emx, sep="\t", index_col=0)
	netlist = pd.read_csv(args.netlist, sep="\t")

	print("Loaded expression matrix (%d genes, %d samples)" % emx.shape)
	print("Loaded netlist (%d edges)" % len(netlist.index))

	# setup plot limits
	if args.pw_scale:
		limits = (emx.min().min(), emx.max().max())
	else:
		limits = None

	# plot cluster distribution
	if args.clusdist:
		print("Plotting cluster distribution...")
		plot_clusdist(netlist, output_dir=args.output_dir)

	# plot correlation distribution
	if args.corrdist:
		print("Plotting correlation distribution...")
		plot_corrdist(netlist, output_dir=args.output_dir)

	# plot gene coverage
	if args.coverage:
		print("Plotting gene coverage...")
		plot_coverage(netlist, output_dir=args.output_dir)

	# plot pairwise plots
	if args.pairwise:
		print("Plotting pairwise plots...")
		plot_pairwise(emx, netlist, output_dir=args.output_dir, limits=limits, range_args=args.pw_range)
