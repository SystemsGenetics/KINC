import argparse
import matplotlib.pyplot as plt
import os
import pandas as pd
import scipy.stats
import seaborn as sns



def load_ematrix(filename):
	print("Loading expression matrix...")

	# read expression matrix from file
	emx = pd.read_table(filename, index_col=0)

	# remove NAN columns
	emx = emx.dropna(axis=1, how="all")

	# transpose to make samples row-wise
	emx = emx.T

	print("Loaded expression matrix (%d samples, %d genes)" % emx.shape)

	return emx



def load_netlist(filename):
	print("Loading netlist...")

	netlist = pd.read_table(filename)

	print("Loaded netlist (%d edges)" % len(netlist.index))

	return netlist



if __name__ ==  "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("-e", "--emx", required=True, help="expression matrix file", dest="EMX")
	parser.add_argument("-n", "--netlist", required=True, help="netlist file", dest="NETLIST")
	parser.add_argument("-o", "--output", required=True, help="output directory", dest="OUTPUT")
	parser.add_argument("-s", "--scale", action="store_true", help="use a uniform global scale", dest="SCALE")

	args = parser.parse_args()

	# load input data
	emx = load_ematrix(args.EMX)
	netlist = load_netlist(args.NETLIST)

	# setup plot limits
	if args.SCALE:
		limits = (emx.min(0).min(), emx.max(0).max())
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
		samples = edge["Samples"]

		print(x, y, k)

		# extract pairwise data
		data = emx[[x, y]].dropna(axis=0, how="any")

		# highlight samples in the edge
		colors = []
		for s in samples:
			if s == "9":
				continue
			elif s == "1":
				colors.append("r")
			else:
				colors.append("k")

		# create figure
		plt.subplots(1, 2, sharex=True, sharey=True, figsize=(10, 5))

		# create density plot
		plt.subplot(121)
		plt.xlim(limits)
		plt.ylim(limits)
		sns.kdeplot(data[x], data[y], shade=True, shade_lowest=False)

		# create scatter plot
		plt.subplot(122)
		plt.title("k=%d, samples=%d, r=%0.2f" % (k, edge["Cluster_Samples"], edge["sc"]))
		plt.xlim(limits)
		plt.ylim(limits)
		plt.xlabel(x)
		plt.ylabel(y)
		plt.scatter(data[x], data[y], color="w", edgecolors=colors)

		# save plot to file
		plt.savefig("%s/%s_%s_%d.png" % (args.OUTPUT, x, y, k))
		plt.close()
