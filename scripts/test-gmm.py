import matplotlib.pyplot as plt
import pandas as pd
import sklearn.mixture
import sys



if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "usage: python test-gmm.py [infile]"
		sys.exit(1)

	# load data
	emx = pd.read_csv(sys.argv[1], sep="\t")

	# iterate through each pair
	for i in xrange(len(emx.index)):
		for j in xrange(i):
			# extract pairwise data
			X = emx.iloc[[i, j]].dropna(axis=1, how="any")
			X = X.values.T
			N = X.shape[0]

			# make sure there are enough samples
			min_K = 0
			min_crit = float("inf")

			if N >= 30:
				# initialize clustering models
				models = [sklearn.mixture.GaussianMixture(n_components=n+1) for n in xrange(5)]

				# identify number of clusters
				for k, model in enumerate(models):
					# fit model
					model.fit(X)

					# save the best model
					crit = model.aic(X)
					if crit < min_crit:
						min_K = len(model.weights_)
						min_crit = crit

				# plot clustering results
				plt.subplots(1, len(models), True, True, figsize=(5 * len(models), 5))

				for k, model in enumerate(models):
					K = len(model.weights_)
					crit = model.aic(X)
					y = model.predict(X)

					plt.subplot(1, len(models), k + 1)
					plt.scatter(X[:, 0], X[:, 1], s=20, c=y, cmap="brg")
					plt.title("N = %d, K = %d, crit = %g" % (N, K, crit))
					plt.xlabel(emx.index[i])
					plt.ylabel(emx.index[j])

				plt.show()

			print "%4d %4d: N = %4d, K = %4d" % (i, j, N, min_K)
