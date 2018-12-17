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

	return X, y



def compute_gmm(X, n_components):
	# initialize clustering model
	model = sklearn.mixture.GaussianMixture(n_components)

	# fit clustering model
	model.fit(X)

	# save clustering results
	K = n_components
	y = model.predict(X)

	# compute criterion value
	crit = model.bic(X)

	# print results
	print("%4d %4d: %8s: K = %4d, crit = %g" % (i, j, "GMM", K, crit))

	return K, y, crit



def compute_vbgmm(X, n_components, weight_concentration_prior, weight_threshold):
	# initialize clustering model
	model = sklearn.mixture.BayesianGaussianMixture(n_components, weight_concentration_prior=weight_concentration_prior)

	# fit clustering model
	model.fit(X)

	print("".join(["%8.3f" % w for w in model.weights_]))

	# compute number of effective components
	K = sum([(w > weight_threshold) for w in model.weights_])

	# save clustering results
	y = model.predict(X)

	# print results
	print("%4d %4d: %8s: y_0 = %g, K = %d" % (i, j, "VBGMM", weight_concentration_prior, K))

	return K, y



if __name__ == "__main__":
	if len(sys.argv) != 2:
		print("usage: python test-vbgmm.py [infile]")
		sys.exit(1)

	# define parameters
	min_expression = float("-inf")
	min_samples = 30
	min_clusters = 1
	max_clusters = 5
	weight_concentration_priors = [1e-6, 1e-3, 1e0, 1e3, 1e6]
	weight_threshold = 0.05

	# load data
	emx = pd.read_table(sys.argv[1])

	# iterate through each pair
	for i in range(len(emx.index)):
		for j in range(i):
			# extract clean pairwise data
			X, y = fetch_pair(emx, i, j, min_expression)
			X = X[y == 0]

			if len(X) < min_samples:
				continue

			# compute Gaussian mixture models
			gmms = []

			for n_components in range(min_clusters, max_clusters + 1):
				gmms.append(compute_gmm(X, n_components))

			# compute variational Bayesian Gaussian mixture models
			vbgmms = []

			for weight_concentration_prior in weight_concentration_priors:
				vbgmms.append(compute_vbgmm(X, max_clusters, weight_concentration_prior, weight_threshold))

			# plot comparison of GMMs and VBGMMs
			rows, cols = 2, max(len(gmms), len(vbgmms))
			plt.figure(figsize=(5 * cols, 5 * rows))

			for k in range(len(gmms)):
				K, y, crit = gmms[k]

				plt.subplot(rows, cols, k + 1)
				plt.scatter(X[:, 0], X[:, 1], s=20, c=y, cmap="brg")
				plt.title("GMM: K = %d, crit = %g" % (K, crit))
				plt.xlabel(emx.index[i])
				plt.ylabel(emx.index[j])

			for k in range(len(vbgmms)):
				K, y = vbgmms[k]

				plt.subplot(rows, cols, cols + k + 1)
				plt.scatter(X[:, 0], X[:, 1], s=20, c=y, cmap="brg")
				plt.title("VBGMM: y_0 = %.0e, K = %d" % (weight_concentration_priors[k], K))
				plt.xlabel(emx.index[i])
				plt.ylabel(emx.index[j])

			plt.show()
