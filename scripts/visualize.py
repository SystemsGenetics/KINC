import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy.stats
import seaborn as sns



if __name__ ==  "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("-e", "--emx", required=True, help="expression matrix file", dest="EMX")
	parser.add_argument("-n", "--netlist", required=True, help="netlist file", dest="NETLIST")
	parser.add_argument("-o", "--output", required=True, help="output directory", dest="OUTPUT")
	parser.add_argument("-s", "--scale", action="store_true", help="use a uniform global scale", dest="SCALE")

	args = parser.parse_args()

	# load input data
	emx = pd.read_table(args.EMX, index_col=0)
	netlist = pd.read_table(args.NETLIST)

	print("Loaded expression matrix (%d genes, %d samples)" % emx.shape)
	print("Loaded netlist (%d edges)" % len(netlist.index))

	# setup plot limits
	if args.SCALE:
		limits = (emx.min().min(), emx.max().max())
	else:
		limits = None

	# initialize output directory
	if not os.path.exists(args.OUTPUT):
		os.mkdir(args.OUTPUT)

	# iterate through each network edge
	for idx in netlist.index:
		edge = netlist.iloc[idx]
		x = edge["Source"]
		y = edge["Target"]
		k = edge["Cluster"]

		print(x, y, k)

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
		plt.xlim(limits)
		plt.ylim(limits)
		sns.kdeplot(data[0], data[1], shade=True, shade_lowest=False)

		# create scatter plot
		plt.subplot(122)
		plt.title("k=%d, samples=%d, spearmanr=%0.2f" % (k, edge["Cluster_Samples"], r))
		plt.xlim(limits)
		plt.ylim(limits)
		plt.xlabel(x)
		plt.ylabel(y)
		plt.scatter(data[0], data[1], color="w", edgecolors=colors)

		# save plot to file
		plt.savefig("%s/%s_%s_%d.png" % (args.OUTPUT, x, y, k))
		plt.close()
