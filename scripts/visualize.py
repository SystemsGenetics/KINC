import argparse
import matplotlib.pyplot as plt
import os
import pandas as pd
import scipy.stats
import seaborn as sns



def load_ematrix(filename):
	print "Loading expression matrix..."

	# read expression matrix from file
	emx = pd.read_csv(filename, sep="\t", index_col=0)

	# remove NAN columns
	emx = emx.dropna(axis=1, how="all")

	# transpose to make samples row-wise
	emx = emx.T

	print "Loaded expression matrix (%d samples, %d genes)" % emx.shape

	return emx



def load_netlist(filename):
	print "Loading netlist..."

	netlist = pd.read_csv(filename, sep="\t")

	print "Loaded netlist (%d edges)" % len(netlist.index)

	return netlist



if __name__ ==  "__main__":
	# define plotting methods
	METHODS = {
		"kde": sns.kdeplot,
		"scatter": plt.scatter
	}

	# parse command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("-e", "--emx", required=True, help="expression matrix file", dest="EMX")
	parser.add_argument("-m", "--method", default="scatter", choices=["kde", "scatter"], help="plotting method", dest="METHOD")
	parser.add_argument("-n", "--netlist", required=True, help="netlist file", dest="NETLIST")
	parser.add_argument("-o", "--output", required=True, help="output directory", dest="OUTPUT")
	parser.add_argument("-s", "--scale", action="store_true", help="use a uniform global scale", dest="SCALE")

	args = parser.parse_args()

	# load input data
	emx = load_ematrix(args.EMX)
	netlist = load_netlist(args.NETLIST)

	# select plotting method
	joint_func = METHODS[args.METHOD]

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

		print x, y, k

		# highlight samples in the edge
		colors = []
		for s in samples:
			if s == "9":
				continue
			elif s == "1":
				colors.append("r")
			else:
				colors.append("k")

		# create joint plot
		g = sns.JointGrid(x=x, y=y, data=emx, xlim=limits, ylim=limits)
		g = g.plot_joint(joint_func, color="w", edgecolor=colors)
		g = g.plot_marginals(sns.distplot, kde=False, color="k")

		# add correlation (might fail)
		try:
			g = g.annotate(scipy.stats.pearsonr)
		except:
			pass

		# save plot to file
		g.savefig("%s/%s_%s_%d.png" % (args.OUTPUT, x, y, k))
		plt.close()
