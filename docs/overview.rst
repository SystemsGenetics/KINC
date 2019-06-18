Overview
========

.. figure:: images/kinc.png
   :alt: KINC Logo

.. image:: https://zenodo.org/badge/71836133.svg
   :target: https://zenodo.org/badge/latestdoi/71836133

The Knowledge Independent Network Construction (KINC) software is a C++ application for constructing a `gene co-expression network (GCN) <https://en.wikipedia.org/wiki/Gene_co-expression_network>`_ from gene expression data organized in a text file called a Gene Expression Matrix (GEM).

How does KINC work?
-------------------

A GEM contain an *n*x*m* matrix of data where *n* is the number of genes and *m* is the number of measurements (or samples) and data is typically normalized and log transformed gene expression data from either RNA-Seq or microarray data. The tool `GEMmaker <https://github.com/SystemsGenetics/GEMmaker>`_ can help create GEMs from RNA-Seq data using tools such as `Hisat2 <https://ccb.jhu.edu/software/hisat2/index.shtml>`_, `Kallisto <https://pachterlab.github.io/kallisto/>`_ or `Salmon <https://combine-lab.github.io/salmon/>`_.

Construction of GCNs traditionally involves several steps that includes, pairwise correlation analysis,  thresholding and module discovery.  Currently, the most popular GCN construction tool, `WGCNA <https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/>`_, performs all of these steps. Although for WGCNA, thresholding and module discovery are a single process. KINC builds from such tools, but uses a new approach to reduce noise that leads to false positives, reduces false negatives, and yields context-specific subgraphs. Often, gene interactions that are specific to a particular experimental condition, tissue type, developmental stage, time series or genotype are missing from the network if they do not represent a major component of variation across all samples. KINC helps tease these out, and creates context-specific subgraph, specific to a given condition. In summary KINC provides the following methods:

**Pairwise Clustering**

- *Gaussian Mixture Models (GMMs)* is performed prior to correlation analysis, modes of gene expression are identified. These modes potentially represent condition-specific gene relationships and become `edges` in teh network if significant.

**Pairwise Correlation Analysis**

In KINC, GMM clusters undergo correlation analysis independent of one another. This is different from traditional approaches where all samples, which may not meet the assumptions of the correlation test, are analyzed together.  Methods include:

- *Pearson*
- *Spearman*

.. note::

  KINC does not implement Mutual Information, another common "association" method for GCN construction as literature has shown no significant improvement over methods such as Pearson's or Spearman.

**Thresholding**

- *Correlation Power Analysis* removes non-significant edges. Correlations from GMM clusters that had insufficient samples to "trust" a given correlation are removed.
- *Power-law* thresholding iterates through decreasing correlation scores, and checks if the resulting network follows a power-law distribution based on node degree. A network that follows a power-law distribution is known as a `scale-free` network. A threshold can be identified at a given correlation score where the resulting network does not appear scale-free.
- *Random Matrix Theory (RMT)* also iterates through decreasing correlation scores, the resulting similarity matrix is examined for properties of random matricies. If the resulting matrix ceases to look non-random a threshold is identified.

.. note::

  A user of KINC may decide to apply one or more of the above thresholding approaches.

**Module Discovery**

How to run KINC?
----------------

KINC can be run on a stand-alone Linux desktop or a High Performance Computing (HPC) cluster.  Traditional network construction (without GMMs) can easily be performed on a stand-alone machine.  However, use of GMMs requires GPUs, and as the size of the GEM grows larger, KINC requires multiple GPUs.  Therefore, it is recommended to use KINC on an HPC system.  Not everyone, unfortunately, has access or is comfortable using HPC systems but due to the computational time, execution on HPC is the only way to create context-specific subgraphs.

KINC provides a graphical interface as well as a command-line.  The GUI is useful for fast running tasks, such as importing of the GEM file, thresholding and network file export.  But the command-line version should be used for GMMs and correlation analysis.  See the :doc:`usage` page for instructions.

How was KINC created?
---------------------

KINC is built with `ACE <https://github.com/SystemsGenetics/ACE>`__, a framework which provides mechanisms for large-scale heterogeneous computing and data management. As such, KINC can be run in a variety of compute configurations, including single-CPU / single-GPU and multi-CPU / multi-GPU, and KINC uses its own binary file formats to represent the data objects that it consumes and produces. Each of these binary formats can be exported to a plain-text format for use in other applications.
