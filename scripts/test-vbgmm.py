import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn.mixture
import sys



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



def compute_gmm(X, y, min_samples, min_clusters, max_clusters, criterion):
	# extract clean pairwise data
	mask = (y == 0)
	X_clean = X[mask]
	y_out = np.copy(y)

	# make sure there are enough samples
	N = X_clean.shape[0]
	K = 0

	if N >= min_samples:
		# initialize clustering models
		models = [sklearn.mixture.GaussianMixture(K) for K in range(min_clusters, max_clusters+1)]
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
				y_out[mask] = model.predict(X_clean)

	print("%4d %4d: %8s: N = %4d, K = %4d" % (i, j, "GMM", N, K))

	return K, y_out



def compute_vbgmm(X, y, min_samples, n_components, weight_concentration_prior, weight_threshold):
	# extract clean pairwise data
	mask = (y == 0)
	X_clean = X[mask]
	y_out = np.copy(y)

	# make sure there are enough samples
	N = X_clean.shape[0]
	K = 0

	if N >= min_samples:
		# initialize clustering model
		model = sklearn.mixture.BayesianGaussianMixture(n_components, weight_concentration_prior=weight_concentration_prior)

		# fit clustering model
		model.fit(X_clean)

		print("".join(["%8.3f" % w for w in model.weights_]))

		# compute number of effective components
		K = sum([(w > weight_threshold) for w in model.weights_])

		# plot clustering results
		y_out[mask] = model.predict(X_clean)

	print("%4d %4d: %8s: N = %4d, K = %4d" % (i, j, "VBGMM", N, K))

	return K, y_out



if __name__ == "__main__":
	if len(sys.argv) != 2:
		print("usage: python test-vbgmm.py [infile]")
		sys.exit(1)

	# define parameters
	min_expression = float("-inf")
	min_samples = 30
	min_clusters = 1
	max_clusters = 5
	criterion = "bic"
	weight_concentration_priors = [1e-3, 1e0, 1e3]
	weight_threshold = 0.05

	# load data
	emx = pd.read_table(sys.argv[1])

	# iterate through each pair
	for i in range(len(emx.index)):
		for j in range(i):
			# fetch pairwise input data
			X, y = fetch_pair(emx, i, j, min_expression)

			models = []

			# compute Gaussian mixture model
			models.append(compute_gmm(X, y, min_samples, min_clusters, max_clusters, criterion))

			# compute variational Bayesian Gaussian mixture model
			for weight_concentration_prior in weight_concentration_priors:
				models.append(compute_vbgmm(X, y, min_samples, max_clusters, weight_concentration_prior, weight_threshold))

			# plot comparison of GMM and VBGMM
			if np.any([K > 0 for K, y in models]):
				plt.figure(figsize=(5 * len(models), 5))

				for i in range(len(models)):
					K, y = models[i]

					if i == 0:
						title = "GMM: K = %d" % (K)
					else:
						title = "VBGMM (y_0 = %g): K = %d" % (weight_concentration_priors[i - 1], K)

					plt.subplot(1, len(models), i + 1)
					plt.scatter(X[:, 0], X[:, 1], s=20, c=y, cmap="brg")
					plt.title(title)
					plt.xlabel(emx.index[i])
					plt.ylabel(emx.index[j])

				plt.show()
