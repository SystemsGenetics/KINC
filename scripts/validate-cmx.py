import argparse
import numpy as np
import pandas as pd



def pairwise_error(pair_true, pair_test, K, column):
	error = 0.0

	for x_true, x_test in zip(pair_true[column], pair_test[column]):
		error += abs(x_true - x_test) / K

	return error



if __name__ ==  "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument(dest="cmx_true", help="true correlation file")
	parser.add_argument(dest="cmx_test", help="test correlation file")

	args = parser.parse_args()

	# load input data
	names = [
		"Source",
		"Target",
		"Cluster_Index",
		"Num_Clusters",
		"Num_Samples",
		"Similarity_Score",
		"Samples"
	]
	cmx_true = pd.read_csv(args.cmx_true, sep="\t", names=names, index_col=[0, 1, 2])
	cmx_test = pd.read_csv(args.cmx_test, sep="\t", names=names, index_col=[0, 1, 2])

	# compore number of pairs
	print("Number of pairs (true): %d" % len(cmx_true.index))
	print("Number of pairs (test): %d" % len(cmx_test.index))

	# get list of all pairs
	pairs_true = [(x, y) for (x, y, k) in cmx_true.index]
	pairs_test = [(x, y) for (x, y, k) in cmx_test.index]
	pairs = list(set(pairs_true + pairs_test))
	pairs.sort()

	# compute error between each pair
	n_shared_edges = 0
	error_K = 0.0
	error_N = 0.0
	error_r = 0.0
	error_S = 0.0

	for x, y in pairs:
		# extract pair from each cmx
		pair_true = cmx_true.loc[(x, y)]
		pair_test = cmx_test.loc[(x, y)]

		# compute error in number of clusters
		K_true = len(pair_true)
		K_test = len(pair_test)

		error_K += abs(K_true - K_test) / len(pairs)

		# report errors
		if K_true != K_test:
			print("%4d %4d: %d != %d" % (x, y, K_true, K_test))

		# use smaller K for cluster-wise comparisons
		K = min(K_true, K_test)

		# update number of shared edges
		n_shared_edges += K

		# compute error in sample size
		error_N += pairwise_error(pair_true, pair_test, K, "Num_Samples") / len(pairs)

		# compute error in similarity score
		error_r += pairwise_error(pair_true, pair_test, K, "Similarity_Score") / len(pairs)

		# compute error in sample mask
		error_S_pair = 0.0

		for S_true, S_test in zip(pair_true["Samples"], pair_test["Samples"]):
			error_S_pair += sum([(s_true != s_test) for s_true, s_test in zip(S_true, S_test)]) / len(S_true) / K

		error_S += error_S_pair / len(pairs)

	print("Number of shared edges: %8d" % (n_shared_edges))
	print("")
	print("Error summary:")
	print("  Number of clusters:   %8.3f" % (error_K))
	print("  Cluster size:         %8.3f" % (error_N))
	print("  Correlation:          %8.3f" % (error_r))
	print("  Sample mask:          %8.3f" % (error_S))
