import argparse
import numpy as np
import pandas as pd



def pairwise_error(pair_true, pair_test, K, column_idx):
	error = 0.0

	for k in range(K):
		x_true = pair_true.iloc[k, column_idx]
		x_test = pair_test.iloc[k, column_idx]

		error += abs(x_true - x_test) / K

	return error



if __name__ ==  "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument(dest="cmx_true", help="true correlation file")
	parser.add_argument(dest="cmx_test", help="test correlation file")

	args = parser.parse_args()

	# load input data
	cmx_true = pd.read_csv(args.cmx_true, sep="\t", header=None, index_col=False)
	cmx_test = pd.read_csv(args.cmx_test, sep="\t", header=None, index_col=False)

	# compore number of pairs
	print("Number of pairs (true): %d" % len(cmx_true.index))
	print("Number of pairs (test): %d" % len(cmx_test.index))

	# get list of all pairs
	pairs_true = [(cmx_true.iloc[idx, 0], cmx_true.iloc[idx, 1]) for idx in cmx_true.index]
	pairs_test = [(cmx_test.iloc[idx, 0], cmx_test.iloc[idx, 1]) for idx in cmx_test.index]
	pairs = list(set(pairs_true + pairs_test))

	pairs.sort()

	# compute pairwise statistics
	n_shared_edges = 0
	error_K = 0.0
	error_N_c = 0.0
	error_N_m = 0.0
	error_N_t = 0.0
	error_N_o1 = 0.0
	error_N_o2 = 0.0
	error_r = 0.0
	error_S = 0.0

	for idx in pairs:
		# extract pair from each cmx
		pair_true = cmx_true.loc[(cmx_true[0] == idx[0]) & (cmx_true[1] == idx[1])]
		pair_test = cmx_test.loc[(cmx_test[0] == idx[0]) & (cmx_test[1] == idx[1])]

		# compute error in number of clusters
		K_true = 0 if pair_true.empty else pair_true.iloc[0, 3]
		K_test = 0 if pair_test.empty else pair_test.iloc[0, 3]

		error_K += abs(K_true - K_test) / len(pairs)

		# report errors
		if K_true != K_test:
			print("%4d %4d: %d != %d" % (idx[0], idx[1], K_true, K_test))

		# use smaller K for cluster-wise comparisons
		K = min(K_true, K_test)

		# update number of shared edges
		n_shared_edges += K

		# compute error in clean sample size
		error_N_c += pairwise_error(pair_true, pair_test, K, 4) / len(pairs)

		# compute error in missing sample size
		error_N_m += pairwise_error(pair_true, pair_test, K, 5) / len(pairs)

		# compute error in thresholded sample size
		error_N_t += pairwise_error(pair_true, pair_test, K, 6) / len(pairs)

		# compute error in thresholded sample size
		error_N_o1 += pairwise_error(pair_true, pair_test, K, 7) / len(pairs)

		# compute error in thresholded sample size
		error_N_o2 += pairwise_error(pair_true, pair_test, K, 8) / len(pairs)

		# compute error in correlation
		error_r += pairwise_error(pair_true, pair_test, K, 9) / len(pairs)

		# compute error in sample mask
		error_S_pair = 0.0

		for k in range(K):
			S_true = pair_true.iloc[k, 10]
			S_test = pair_test.iloc[k, 10]

			error_S_pair += sum([(s_true != s_test) for s_true, s_test in zip(S_true, S_test)]) / len(S_true) / K

		error_S += error_S_pair / len(pairs)

	print("Number of shared edges:     %8d" % (n_shared_edges))
	print("")
	print("Error summary:")
	print("  Number of clusters:       %8.3f" % (error_K))
	print("  Clean sample size:        %8.3f" % (error_N_c))
	print("  Missing sample size:      %8.3f" % (error_N_m))
	print("  Thresholded sample size:  %8.3f" % (error_N_t))
	print("  Pre-outlier sample size:  %8.3f" % (error_N_o1))
	print("  Post-outlier sample size: %8.3f" % (error_N_o2))
	print("  Correlation:              %8.3f" % (error_r))
	print("  Sample mask:              %8.3f" % (error_S))
