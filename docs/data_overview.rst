Data Format Overview
====================

KINC generates several custom "data objects" stored as binary files.  These files represent the data types that KINC consumes and produces. Additionally, KINC can convert these objects to and from more conventional plain-text formats. Each data object and plain-text format is described here.

Gene Expression Matrix (GEM)
----------------------------
The GEM is the starting point for all network construction.  It is a tab-delimited text file containing an `n` x `m` matrix of gene expression values where rows (or lines in the file) are genes, columns are samples and elements are the expression values for a gene in a given sample.

The GEM may have a header line. But the header line must only contain the sample names. Each subsequent line will always start with the gene name, followed by the expression values.  Therefore, the header line, if present, should always have one less element then every other row.

.. note::

  The header line, if present, should always have one less element then every other row. Missing (NAN) values are often represented by a special token such as "NA" or "0.00".

A very simple example is shown below:

.. code::

	Sample1	Sample2	Sample3	Sample4
	Gene1	0.523	0.991	0.421	NA
	Gene2	NA	7.673	3.333	9.103
	Gene3	4.444	5.551	NA	0.013

.. note::

  The GEM is imported into KINC and stored in binary format as an Expression Matrix file (described below).  KINC can also export the expression matrix back into its text-based GEM format.

.. _amx-reference-label:

Sample Annotation Matrix (AMX)
------------------------------
The annotation matrix is a tab-delimited or comma-separated file that contains metadata about the RNA-seq samples used to construct the network.  The first column should list the sample names (the same as in the GEM) and each additional column should contain information about the sample such as experimental condtions (e.g. Treatment, Tissue, Developmental Stage, Sampling Time, Genotype, etc.) which are usually categorical data. The file can also contain phenotype information that may have been collected from the individuals from which the samples were taken.

.. note::

  The annotation matrix is only needed if condition-specific subgraphs are wanted.

Binary Output Files
-------------------

Expression Matrix (EMX)
~~~~~~~~~~~~~~~~~~~~~~~
The Expression Matrix (EMX) file is a data matrix, with genes as rows and samples as columns. It is constructed from the tab-delimited plain-text input Gene Expression Matrix (GEM) file. It contains the matrix data in binary format as well as the gene names and sample names. EMX files have the ``.emx`` extension.

Correlation Matrix (CMX)
~~~~~~~~~~~~~~~~~~~~~~~~
The Correlation Matrix (CMX) file is created during the ``similarity`` step and is used to represent a similarity matrix, which is a pairwise matrix of similarity (e.g. correlation) scores for each gene pair using data from an expression matrix (EMX). For traditional network construction, the similarity matrix will contain a single similarity score. For the GMM approach, the similarity matrix will contain multiple scores, one for each cluster identifed using GMMs. The CMX uses a sparse matrix format--it only stores the gene pairs which have correlation data. Otherwise, the file would grow too large for most file systems. CMX files have the ``.cmx`` extension.

Cluster Composition Matrix (CCM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The Cluster Composition Matrix (CCM) file is created during the ``similarity`` step only when using the GMM approach to network construction.  This maintains a matrix identical to the CMX file but instead of similarity scores for each clsuter, it stores a "sample composition string".  This string is a series of numbers, strung together indicating which samples belong to the cluster.  If a 1 appears in the first position, it indicates the the first sample in the EMX file belongs to the cluster. If a 0 appears it indicates the sample does not.  Several other numbers are present to indicate missing or outlier samples. The same is true for the 2nd, 3rd and up to the nth samples. The following indicates the meaning of each number:

- ``0``: The sample is not present in the cluster.
- ``1``: The sample is present in the cluster.
- ``6``: The sample was removed because the expression level in one gene was below the allowed expression threshold.
- ``7``: The sample was removed as an outlier prior to GMMs.
- ``8``: The sample was removed as an outlier after the GMM cluster was identified.
- ``9``: The sample was ignored because one of the genes was missing expression in the sample.


The CCM also uses a sparse matrix format and has a ``.ccm`` extension.

Condition Specific Matrix (CSM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The Condition Specific Matrix (CSM) file is created when the ``cond-test`` thresholding function is used. This thresholding approach can only be used if the GMM approach was used. It requires an “Annotation Matrix” as input, which is a tab-delimited file where the rows are samples and each column contains some feature about the sample, such as experimental conditions or phenotypic values. The file will contain a matrix identical to the CMX matrix but instead of similarity scores, it stores p-values for the association of selected feature from the annotation matrix. The CSM also uses a sparse matrix format and has a ``.csm`` extension.

.. _plain-text-reference-label:

Plain-text Output Files
-----------------------

Network File
~~~~~~~~~~~~
There are two types of plan-text network files. First is a tab-delimited file where each line contains information about a single edge in the network. The file contains the following columns:

- ``Source``:  The first gene in the edge
- ``Target``:  Thes second gene in the edge
- ``Similarity_Score``:  The similarity score value.
- ``Interaction``: The interaction type. This always is ``co`` for co-expression.
- ``Cluster_Index``: The unique cluster index (starting from zero)
- ``Cluster_Size``: The total number of samples in the cluster.
- ``Samples``:  The sample composition string.

Additionally, if the ``cond-test`` function was performed, a series of additional columns will be present containing the p-values for each test performed.

The following is a sample line from a network file:

.. code:: bash

	Source	Target  Similarity_Score  Interaction	Cluster_Index	Cluster_Size Samples
	Gene1	Gene2	0.979	co	0	30	1199991911111161111111611161111111111770080000000

Additionally, KINC does support creation of a "minimal" plain-text format, which does not contain the sample string or summary statistics. This format is useful for inspecting large networks quickly. The following is a sample line of a minimal network file:

.. code:: bash

	Source	Target	sc	Cluster	Num_Clusters
	Gene1	Gene2	0.979	0	1

The second major network file format is the GraphML format. This is a common XML format used for representing networks. The following is an example snippet of a GraphML file generated by KINC:

.. code:: XML

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

Correlation Matrix
~~~~~~~~~~~~~~~~~~
A plain-text correlation matrix is a representation of a sparse matrix where each line is a correlation. It includes the pairwise index, correlation value, sample composition string, and several other summary statistics.  The following is a sample line from the correlation matrix file:

.. code:: bash

	0	1	0	1	30	5	2	1	3	0.979	1199991911111161111111611161111111111770080000000
