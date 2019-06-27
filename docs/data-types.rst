Data Types
==========

KINC uses several custom "data objects" to represent the data types that it consumes and produces. Additionally, KINC can convert these objects to and from more conventional plain-text formats. Each data object and plain-text format is described here.

Data Objects
------------

Expression Matrix (EMX)
~~~~~~~~~~~~~~~~~~~~~~~

The Expression Matrix (EMX) object is a data matrix, with genes as rows and samples as columns. It contains the matrix data in binary format as well as the gene names and sample names. EMX files have the ``.emx`` extension.

Cluster Composition Matrix (CCM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Cluster Composition Matrix (CCM) object is one of two files used to represent a similarity matrix, which is a pairwise matrix of similarity scores for each gene pair in an expression matrix. The similarity matrix can contain multiple scores for a single pair, which correspond to clusters of pairwise samples, and the CCM contains the cluster labels, or "sample strings" for each pairwise cluster. The CCM uses a sparse matrix format -- it only stores the gene pairs which have cluster data. CCM files have the ``.ccm`` extension.

Correlation Matrix (CMX)
~~~~~~~~~~~~~~~~~~~~~~~~

The Correlation Matrix (CMX) object is one of the two files used to represent a similarity matrix, which is a pairwise matrix of similarity scores for each gene pair in an expression matrix. The similarity matrix can contain multiple scores for a single pair, which correspond to clusters of pairwise samples, and the CMX contains the correlations for each pairwise cluster. The CMX uses a sparse matrix format -- it only stores the gene pairs which have correlation data. CMX files have the ``.cmx`` extension.

Plain-text Formats
------------------

Expression Matrix
~~~~~~~~~~~~~~~~~

A plain-text expression matrix is a data matrix with each row on a line, each value separated by whitespace, and the first row and column containing the row names and column names, respectively. Missing (NAN) values are represented by a special token such as "NA". Expression matrix files have the ``.emx.txt`` extension.

.. code:: bash

		Sample1	Sample2	Sample3	Sample4
	Gene1	0.523	0.991	0.421	NA
	Gene2	NA	7.673	3.333	9.103
	Gene3	4.444	5.551	NA	0.013

Correlation Matrix
~~~~~~~~~~~~~~~~~~

A plain-text correlation matrix is a sparse matrix format where each line is a correlation that includes the pairwise index, correlation value, sample mask, and several other summary statistics.

.. code:: bash

	0	1	0	1	30	5	2	1	3	0.979	1199991911111161111111611161111111111770080000000

Co-expression Network
~~~~~~~~~~~~~~~~~~~~~

A co-expression network is an edge list where each line is an edge between two genes. It is very similar to the plain-text correlation matrix.

.. code:: bash

	Source	Target	sc	Interaction	Cluster	Num_Clusters	Cluster_Samples	Missing_Samples	Cluster_Outliers	Pair_Outliers	Too_Low	Samples
	Gene1	Gene2	0.979	co	0	1	30	5	2	1	3	1199991911111161111111611161111111111770080000000

There is also a "minimal" plain-text format, which does not contain the sample string or summary statistics. This format is useful for inspecting large networks quickly.

.. code:: bash

	Source	Target	sc	Cluster	Num_Clusters
	Gene1	Gene2	0.979	0	1

There is also the GraphML format, which is an XML format for representing networks.

.. code:: bash

	<?xml version="1.0" encoding="UTF-8"?>
	<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
	         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
	         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
		<graph id="G" edgedefault="undirected">
			<node id="Gene1"/>
			<node id="Gene2"/>
			<edge source="Gene1" target="Gene2" samples="1199991911111161111111611161111111111770080000000"/>
		</graph>
	</graphml>
