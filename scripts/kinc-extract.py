import argparse
import pandas as pd



if __name__ == "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("--emx", help="expression matrix file", required=True)
	parser.add_argument("--cmx", help="correlation matrix file", required=True)
	parser.add_argument("--output", help="output net file", required=True)
	parser.add_argument("--mincorr", help="minimum absolute correlation threshold", type=float, default=0)
	parser.add_argument("--maxcorr", help="maximum absolute correlation threshold", type=float, default=1)

	args = parser.parse_args()

	# load data
	emx = pd.read_csv(args.emx, sep="\t")
	cmx = pd.read_table(args.cmx, header=None, names=[
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
	cmx = cmx[(args.mincorr <= abs(cmx["sc"])) & (abs(cmx["sc"]) <= args.maxcorr)]

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
	cmx.to_csv(args.output, sep="\t", index=False)
