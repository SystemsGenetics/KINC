import argparse
import pandas as pd
import sklearn.datasets



if __name__ == "__main__":
	# parse command-line arguments
	parser = argparse.ArgumentParser(description="Create a synthetic gene expression dataset")
	parser.add_argument("--n-samples", help="number of samples", type=int, default=100)
	parser.add_argument("--n-genes", help="number of genes", type=int, default=20)
	parser.add_argument("--n-classes", help="number of classes", type=int, default=2)
	parser.add_argument("--dataset", help="name of dataset file", default="GEM.txt")
	parser.add_argument("--transpose", help="transpose the dataset to be (genes x samples)", action="store_true")

	args = parser.parse_args()

	# create synthetic dataset
	n_informative = args.n_genes // 10
	n_redundant = args.n_genes - n_informative

	X, y = sklearn.datasets.make_classification(args.n_samples, args.n_genes, n_informative=n_informative, n_redundant=n_redundant, n_classes=args.n_classes)

	# initialize gene names, sample names
	X_samples = ["sample-%08d" % i for i in range(args.n_samples)]
	X_genes = ["gene-%06d" % i for i in range(args.n_genes)]

	# initialize dataframe
	X = pd.DataFrame(X, index=X_samples, columns=X_genes)

	# transpose dataset if specified
	if args.transpose:
		X = X.T

	# save dataset to file
	X.to_csv(args.dataset, sep="\t", na_rep="NA", float_format="%.8f")
