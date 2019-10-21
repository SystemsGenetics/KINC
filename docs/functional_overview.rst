Functional Overview
===================

KINC provides two executables: ``kinc`` (the command-line version) and ``qkinc`` (the graphical user interface (GUI) version). Both command-line and GUI version can perform all of the same functionality, except that the command-line version can use the `Message Passing Interface <https://www.open-mpi.org/>`_ (MPI) for inter-process communication on an HPC system. The GUI cannot.

KINC can construct GCNs using a traditional approach (see the :doc:`about` section) or using the new GMM approach.  To this end, a variety of functions for each step in the GCN construction workflow are provided. The names and descriptions of these functions are provided below in the order that they might be used for network construction.

.. note::

  The command-line and GUI versions of KINC provide the same set of functions.  The names of the command-line tool are provided below but the descriptions are the same for both.

Step 1: GEM Import
  - ``import-emx``: Import a gene expression matrix, or GEM, into a binary format. The binary format is readable by other KINC functions.


Step 2: Simlarity Matrix Construction
  - ```similarity``: This function is responsible for creating the similarity matrix and supports both traditional and GMM approaches for GCN construction. It uses the GEM file imported using the ``import-emx`` function.

Step 3: Filtering (optional)
  - ``corrpower``: This function performs `power analysis <https://www.statmethods.net/stats/power.html>`_ and removes edges from the network that have insufficient power. This function is only necessary when the number of samples in a cluster are allowed to be small.  The minimum size cluster can be set using the ``similarity`` step. For small clusters, the Type I and Type II error rates are higher and this function ensures that the only edges that remain in the final network are those with a sufficient alpha (default 0.001) and beta (default 0.8) values.

Step 4: Thresholding
  - ``rmt``: This function provides the Random Matrix Theory (RMT) approach to thresholding a network. It iterates through decreasing correlation scores, the resulting similarity matrix is examined for properties of random matrices. If the resulting matrix ceases to look non-random a threshold is identified. It is not fully compatible with the GMM approach.
  - ``powerlaw``: This thresholds the network at the point that it appears `scale-free`. Most real networks are scale free and a network that follows a power-law distribution is known as a scale-free.  This approach iterates through decreasing correlation scores, and checks if the resulting network follows a power-law distribution based on node degree. A threshold can be identified at a given correlation score where the resulting network does not appear scale-free.
  - ``cond-test``: This function performs a sequence of binomial and regression tests that test if an edge in the network (i.e. cluster identified by GMMs) is condition-specific.  It requires an "Annotation Matrix": a tab-delimited file where the rows are samples and each column contains some feature about the sample, such as experimental conditions or phenotypic values.  This function is **only applicable to the GMM approach** of network construction and is the primary method for identification of condition specific subgraphs.

Step 5: Network Creation
  - ``extract``: Once the network construction steps are completed, this program extracts the network into a tab-delimeted text file or GraphML file for use by other tools such as `Cytoscape <https://cytoscape.org/>`_.

Other Useful Functions
  - ``export-emx``: Export a gene expression matrix, or GEM, from the binary format used by KINC into a tab-delimited format compatible with other tools.
  - ``export-cmx``: Export a correlation matrix (or similarity matrix) from the binary format that KINC uses into a text-based version.
  - ``import-cmx``: Import a correlation matrix (or similarity matrix) that was previously exported by KINC.
