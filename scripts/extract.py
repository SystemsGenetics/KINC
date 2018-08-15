import argparse
import pandas as pd



if __name__ == "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--emx", required=True, help="expression matrix file", dest="EMX")
	parser.add_argument("--cmx", required=True, help="correlation matrix file", dest="CMX")
	parser.add_argument("-o", "--output", required=True, help="output net file", dest="OUTPUT")
	parser.add_argument("--mincorr", type=float, default=0, help="minimum absolute correlation threshold", dest="MINCORR")
	parser.add_argument("--maxcorr", type=float, default=1, help="maximum absolute correlation threshold", dest="MAXCORR")

	args = parser.parse_args()

	# load data
	emx = pd.read_csv(args.EMX, sep="\t")
	cmx = pd.read_csv(args.CMX, sep="\t", header=None, names=[
		"x",
		"y",
		"Cluster",
		"Num_Clusters",
		"Cluster_Samples",
		"Missing_Samples",
		"Cluster_Outliers",
		"Pair_Outliers",
		"Too_Low",
		"sc",
		"Samples"
	])

	# extract correlations within thresholds
	cmx = cmx[(args.MINCORR <= abs(cmx["sc"])) & (abs(cmx["sc"]) <= args.MAXCORR)]

	# insert additional columns used in netlist format
	cmx.insert(len(cmx.columns), "Source", [emx.index[x] for x in cmx["x"]])
	cmx.insert(len(cmx.columns), "Target", [emx.index[y] for y in cmx["y"]])
	cmx.insert(len(cmx.columns), "Interaction", ["co" for idx in cmx.index])

	# reorder columns to netlist format
	cmx = cmx[[
		"Source",
		"Target",
		"sc",
		"Interaction",
		"Cluster",
		"Num_Clusters",
		"Cluster_Samples",
		"Missing_Samples",
		"Cluster_Outliers",
		"Pair_Outliers",
		"Too_Low",
		"Samples"
	]]

	# save output data
	cmx.to_csv(args.OUTPUT, sep="\t", index=False)
