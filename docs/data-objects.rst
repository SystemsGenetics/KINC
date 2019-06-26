Data Objects
============

KINC uses several custom "data objects" to represent the data types that it consumes and produces. Each data object is described here.

Expression Matrix (EMX)
-----------------------

The Expression Matrix (EMX) object is a data matrix, with genes as rows and samples as columns. It contains the matrix data in binary format as well as the gene names and sample names. EMX files have the ``emx`` extension.

Cluster Composition Matrix (CCM)
--------------------------------

The Cluster Composition Matrix (CCM) object is one of two files used to represent a similarity matrix, which is a pairwise matrix of similarity scores for each gene pair in an expression matrix. The similarity matrix can contain multiple scores for a single pair, which correspond to clusters of pairwise samples, and the CCM contains the cluster labels, or "sample strings" for each pairwise cluster. The CCM uses a sparse matrix format -- it only stores the gene pairs which have cluster data. CCM files have the ``ccm`` extension.

Correlation Matrix (CMX)
------------------------

The Correlation Matrix (CMX) object is one of the two files used to represent a similarity matrix, which is a pairwise matrix of similarity scores for each gene pair in an expression matrix. The similarity matrix can contain multiple scores for a single pair, which correspond to clusters of pairwise samples, and the CMX contains the correlations for each pairwise cluster. The CMX uses a sparse matrix format -- it only stores the gene pairs which have correlation data. CMX files have the ``cmx`` extension.
