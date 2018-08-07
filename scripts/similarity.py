import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
import seaborn as sns
import sklearn.cluster
import sklearn.mixture



def create_gmm(n_clusters):
	return sklearn.mixture.GaussianMixture(n_components=n_clusters)



def create_kmeans(n_clusters):
	return sklearn.cluster.KMeans(n_clusters=n_clusters, n_jobs=-1)



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
	# TODO: min expression
	parser.add_argument("--minsamp", type=int, default=30, help="minimum sample size", dest="MINSAMP")
	parser.add_argument("--minclus", type=int, default=1, help="minimum clusters", dest="MINCLUS")
	parser.add_argument("--maxclus", type=int, default=5, help="maximum clusters", dest="MAXCLUS")
	# TODO: criterion type
	# TODO: remove pre-outliers
	# TODO: remove post-outliers
	# TODO: min correlation
	# TODO: max correlation
	parser.add_argument("--visualize", action="store_true", help="whether to visualize results", dest="VISUALIZE")

	args = parser.parse_args()

	# print arguments
	for key, value in vars(args).iteritems():
		print("%s: %s" % (key, str(value)))

	# load data
	emx = pd.read_csv(args.INPUT, sep="\t")
	cmx = open(args.OUTPUT, "w");

	# initialize methods
	clus_method = CLUSTERING_METHODS[args.CLUSMETHOD]
	corr_method = CORRELATION_METHODS[args.CORRMETHOD]

	# iterate through each pair
	for i in xrange(len(emx.index)):
		for j in xrange(i):
			# extract pairwise data
			X = emx.iloc[[i, j]].values.T
			mask = ~np.isnan(X[:, 0]) & ~np.isnan(X[:, 1])

			# initialize labels
			y = np.zeros((X.shape[0],), dtype=int)
			y[~mask] = -9

			# extract clean pairwise data
			X_clean = X[mask]
			N = X_clean.shape[0]

			# peform clustering
			K = 1

			if args.CLUSMETHOD != "none":
				# make sure there are enough samples
				if N >= args.MINSAMP:
					# initialize clustering models
					models = [clus_method(n) for n in xrange(args.MINCLUS, args.MAXCLUS+1)]
					min_crit = float("inf")

					# identify number of clusters
					for model in models:
						# fit model
						model.fit(X_clean)

						# save the best model
						crit = model.bic(X_clean)
						if crit < min_crit:
							min_crit = crit
							K = len(model.weights_)
							y[mask] = model.predict(X_clean)

			print("%4d %4d %d" % (i, j, K))

			# perform correlation
			for k in xrange(K):
				# extract samples in cluster k
				X_k = X[y == k]

				# make sure there are enough samples
				if X_k.shape[0] >= args.MINSAMP:
					# compute correlation
					corr, p = corr_method(X_k[:, 0], X_k[:, 1])

					# compute sample mask
					y_k = np.copy(y)
					y_k[(y_k >= 0) & (y_k != k)] == 0
					y_k[y_k == k] == 1
					y_k[y_k < 0] *= -1

					cmx.write("%d\t%d\t%d\t%d\t%g\t%g\t%s\n" % (i, j, k, K, corr, p, "".join([str(y_i) for y_i in y_k])))

					# plot results
					if args.VISUALIZE:
						sns.jointplot(x=X_k[:, 0], y=X_k[:, 1], kind="reg", stat_func=corr_method)
						plt.show()
