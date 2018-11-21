import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pprint
import scipy.stats
import seaborn as sns
import sklearn.cluster
import sklearn.mixture



def create_gmm(n_clusters):
	return sklearn.mixture.GaussianMixture(n_components=n_clusters)



def create_kmeans(n_clusters):
	return sklearn.cluster.KMeans(n_clusters=n_clusters, n_jobs=-1)



def fetch_pair(emx, i, j, min_expression):
	# extract pairwise data
	X = emx.iloc[[i, j]].values.T

	# initialize labels
	y = np.zeros((X.shape[0],), dtype=int)

	# mark thresholded samples
	y[(X[:, 0] < min_expression) | (X[:, 1] < min_expression)] = -6

	# mark nan samples
	y[np.isnan(X[:, 0]) | np.isnan(X[:, 1])] = -9

	return (X, y)



def mark_outliers(X, labels, k, marker):
	# extract samples in cluster k
	mask = (labels == k)
	x = np.copy(X[mask, 0])
	y = np.copy(X[mask, 1])

	# make sure cluster is not empty
	if len(x) == 0 or len(y) == 0:
		return

	# sort arrays
	x.sort()
	y.sort()

	# compute quartiles and thresholds for each axis
	n = len(x)

	Q1_x = x[n * 1 // 4]
	Q3_x = x[n * 3 // 4]
	T_x_min = Q1_x - 1.5 * (Q3_x - Q1_x)
	T_x_max = Q3_x + 1.5 * (Q3_x - Q1_x)

	Q1_y = y[n * 1 // 4]
	Q3_y = y[n * 3 // 4]
	T_y_min = Q1_y - 1.5 * (Q3_y - Q1_y)
	T_y_max = Q3_y + 1.5 * (Q3_y - Q1_y)

	# mark outliers
	for i in range(len(labels)):
		if labels[i] == k:
			outlier_x = (X[i, 0] < T_x_min or T_x_max < X[i, 0])
			outlier_y = (X[i, 1] < T_y_min or T_y_max < X[i, 1])

			if outlier_x or outlier_y:
				labels[i] = marker



def compute_clustering(X, y, create_model, min_samples, min_clusters, max_clusters, criterion):
	# extract clean pairwise data
	mask = (y == 0)
	X_clean = X[mask]
	N = X_clean.shape[0]

	# make sure there are enough samples
	K = 0

	if N >= min_samples:
		# initialize clustering models
		models = [create_model(K) for K in range(min_clusters, max_clusters+1)]
		min_crit = float("inf")

		# identify number of clusters
		for model in models:
			# fit model
			model.fit(X_clean)

			# compute criterion value
			if criterion == "aic":
				crit = model.aic(X_clean)
			elif criterion == "bic":
				crit = model.bic(X_clean)

			# save the best model
			if crit < min_crit:
				min_crit = crit
				K = len(model.weights_)
				y[mask] = model.predict(X_clean)

	return K, y



def compute_correlation(X, y, k, method, min_samples, visualize):
	# extract samples in cluster k
	X_k = X[y == k]

	# make sure there are enough samples
	if X_k.shape[0] < min_samples:
		return None, None

	# compute correlation
	corr, p = method(X_k[:, 0], X_k[:, 1])

	# plot results
	if visualize:
		sns.jointplot(x=X_k[:, 0], y=X_k[:, 1], kind="reg", stat_func=method)
		plt.show()

	return corr, p



if __name__ == "__main__":
	# define clustering methods
	CLUSTERING_METHODS = {
		"none": None,
		"gmm": create_gmm,
		"kmeans": create_kmeans
	}

	# define correlation methods
	CORRELATION_METHODS = {
		"kendall": scipy.stats.kendalltau,
		"pearson": scipy.stats.pearsonr,
		"spearman": scipy.stats.spearmanr
	}

	# parse command-line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", required=True, help="expression matrix file", dest="INPUT")
	parser.add_argument("-o", "--output", required=True, help="correlation file", dest="OUTPUT")
	parser.add_argument("--clusmethod", default="none", choices=["none", "gmm", "kmeans"], help="clustering method", dest="CLUSMETHOD")
	parser.add_argument("--corrmethod", default="pearson", choices=["kendall", "pearson", "spearman"], help="correlation method", dest="CORRMETHOD")
	parser.add_argument("--minexpr", type=float, default=-float("inf"), help="minimum expression threshold", dest="MINEXPR")
	parser.add_argument("--minsamp", type=int, default=30, help="minimum sample size", dest="MINSAMP")
	parser.add_argument("--minclus", type=int, default=1, help="minimum clusters", dest="MINCLUS")
	parser.add_argument("--maxclus", type=int, default=5, help="maximum clusters", dest="MAXCLUS")
	parser.add_argument("--crit", default="bic", choices=["aic", "bic"], help="model selection criterion", dest="CRITERION")
	parser.add_argument("--preout", action="store_true", help="whether to remove pre-clustering outliers", dest="PREOUT")
	parser.add_argument("--postout", action="store_true", help="whether to remove post-clustering outliers", dest="POSTOUT")
	parser.add_argument("--mincorr", type=float, default=0, help="minimum absolute correlation threshold", dest="MINCORR")
	parser.add_argument("--maxcorr", type=float, default=1, help="maximum absolute correlation threshold", dest="MAXCORR")
	parser.add_argument("--pvalue", type=float, default=float("inf"), help="maximum p-value threshold for correlations", dest="MAXPVALUE")
	parser.add_argument("--visualize", action="store_true", help="whether to visualize results", dest="VISUALIZE")

	args = parser.parse_args()

	# print arguments
	pprint.pprint(vars(args))

	# load data
	emx = pd.read_table(args.INPUT)
	cmx = open(args.OUTPUT, "w");

	# iterate through each pair
	for i in range(len(emx.index)):
		for j in range(i):
			# fetch pairwise input data
			X, y = fetch_pair(emx, i, j, args.MINEXPR)

			# remove pre-clustering outliers
			if args.PREOUT:
				mark_outliers(X, y, 0, -7)

			# perform clustering
			K = 1

			if args.CLUSMETHOD != "none":
				K, y = compute_clustering(X, y, CLUSTERING_METHODS[args.CLUSMETHOD], args.MINSAMP, args.MINCLUS, args.MAXCLUS, args.CRITERION)

			print("%4d %4d %d" % (i, j, K))

			# remove post-clustering outliers
			if K > 1 and args.POSTOUT:
				for k in range(K):
					mark_outliers(X, y, k, -8)

			# perform correlation
			for k in range(K):
				corr, p = compute_correlation(X, y, k, CORRELATION_METHODS[args.CORRMETHOD], args.MINSAMP, args.VISUALIZE)

				# make sure correlation, p-value meets thresholds
				if corr != None and args.MINCORR <= abs(corr) and abs(corr) <= args.MAXCORR and p <= args.MAXPVALUE:
					# compute sample mask
					y_k = np.copy(y)
					y_k[(y_k >= 0) & (y_k != k)] = 0
					y_k[y_k == k] = 1
					y_k[y_k < 0] *= -1

					sample_mask = "".join([str(y_i) for y_i in y_k])

					# compute summary statistics
					num_samples = sum(y_k == 1)
					num_threshold = sum(y_k == 6)
					num_preout = sum(y_k == 7)
					num_postout = sum(y_k == 8)
					num_missing = sum(y_k == 9)

					cmx.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%0.8f\t%s\n" % (i, j, k, K, num_samples, num_missing, num_postout, num_preout, num_threshold, corr, sample_mask))
