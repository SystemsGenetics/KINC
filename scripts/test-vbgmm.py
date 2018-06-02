import matplotlib.pyplot as plt
import pandas as pd
import sklearn.mixture
import sys



if __name__ == "__main__":
	if len(sys.argv) != 2:
		print "usage: python test-vbgmm.py [infile]"
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
			K = 0

			if N >= 30:
				# initialize clustering model
				model = sklearn.mixture.BayesianGaussianMixture(n_components=5, weight_concentration_prior=1e3)

				# fit clustering model
				model.fit(X)

				print "".join(["%8.3f" % w for w in model.weights_])

				# compute number of effective components
				K = sum([1 for w in model.weights_ if w > 0.05])

				# plot clustering results
				y = model.predict(X)
				plt.scatter(X[:, 0], X[:, 1], s=20, c=y, cmap="brg")
				plt.title("N = %d, K = %d" % (N, K))
				plt.xlabel(emx.index[i])
				plt.ylabel(emx.index[j])
				plt.show()

			print "%4d %4d: N = %4d, K = %4d" % (i, j, N, K)
