Functional Overview
===================

KINC provides two executables: ``kinc`` (the command-line version) and ``qkinc`` (the graphical user interface (GUI) version). Both command-line and GUI versions can perform all of the same functionality, except that the command-line version can use the `Message Passing Interface <https://www.open-mpi.org/>`_ (MPI) for inter-process communication on an HPC system. The GUI cannot.  Additionally, there a several scripts that accompany KINC that use `KINC.R <https://github.com/SystemsGenetics/KINC.R>`_ (a supplemental R package for KINC) or custom Python code. As KINC is an evolving tool, these scripts provide functionality that has not yet been incorporated into the KINC binary.

KINC can construct GCNs using a traditional approach (see the :doc:`about` section) or using the new GMM approach.  To this end, a variety of functions for each step in the GCN construction workflow are provided. The names and descriptions of the KINC functions or scripts are provided below in the order that they might be used for network construction.

.. note::

  The command-line and GUI versions of KINC provide the same set of functions.  The names of the command-line tool are provided below but the descriptions are the same for both.

Step 1: GEM Import
  - ``import-emx``: a KINC function that imports a gene expression matrix, or GEM, into a binary format. The binary format is readable by other KINC functions.


Step 2: Simlarity Matrix Construction
  - ``similarity``: a KINC function responsible for creating the similarity matrix and supports both traditional and GMM approaches for GCN construction. It uses the GEM file imported using the ``import-emx`` function.

Step 3: Filtering
  - ``corrpower``: a KINC function that performs `power analysis <https://www.statmethods.net/stats/power.html>`_ and removes edges from the network that have insufficient power. This function is only necessary when the number of samples in a cluster are allowed to be small.  The minimum size cluster can be set using the ``similarity`` step. For small clusters, the Type I and Type II error rates are higher and this function ensures that the only edges that remain in the final network are those with a sufficient alpha (default 0.001) and beta (default 0.8) values.  Skip this function if the minimum cluster size was set to 30 or greater.
  - ``cond-test``: a KINC function that performs a hypergemoetric test to identify an edge in the network (i.e. cluster identified by GMMs) that is condition-specific.  It requires an "Annotation Matrix": a tab-delimited file where the rows are samples and each column contains some feature about the sample, such as experimental conditions or phenotypic values.  Skip this function if GMMs were not used in the ``similarity`` step. This function is **only applicable to the GMM approach** of network construction and is the primary method for identification of condition specific subgraphs.

Step 4: KINC-based Thresholding
  - ``rmt``: This function provides the Random Matrix Theory (RMT) approach to thresholding a network. It iterates through decreasing correlation scores, the resulting similarity matrix is examined for properties of random matrices. If the resulting matrix ceases to look non-random a threshold is identified. It is not fully compatible with the GMM approach. **This is the recommended thresholding approach for traditional network construction.**
  - ``powerlaw``: This thresholds the network at the point that it appears `scale-free`. Most real networks are scale free and a network that follows a power-law distribution is known as a scale-free.  This approach iterates through decreasing correlation scores, and checks if the resulting network follows a power-law distribution based on node degree. A threshold can be identified at a given correlation score where the resulting network does not appear scale-free.

Step 5: Network Extraction
  - ``extract``: Once the network construction steps are completed, this program extracts the network into a tab-delimeted text file or GraphML file for use by other tools such as `Cytoscape <https://cytoscape.org/>`_.

Step 6: Post-Network filtering
  - ``filter-condition-edges.R``:  a KINC.R script that removes two types of biased (or false) edges that can appear in the GMM derived networks:  edges due to lack of differential cluster expression (DCE) or unbalanced missing data between genes in a comparison.  This script is **only applicable to the GMM approach** of network. *Note*: This functionality is in a script until it can be incoporated with the other filter functions of KINC.


Step 7: Post-Network Thresholding
  - ``rank-condition-threshold.R``: a KINC.R script that organizes the edges from a conditional network, ranks them by *p*-value, similarity score, and R:sup:`2` values (for quantitative data) and keeps the top *n* edges.  This is useful for very large condition-specific networks (usually caused by time-series conditions). This script is **only applicable to the GMM approach** of network. *Note*: This functionality is in a script until it can be incoporated with the other threshold functions of KINC.

Step 8:  Network Reports & Visualization
  - ``make-summary-plots.R``:  a KINC.R script that generates publication-quality images describing the relationship between the experimetn conditions, similarity score, *p*-values and R:sup:`2`.  These plots are useful for understanding how condition-specific relationships differ.  This script is **only applicable to the GMM approach** of network.
  - ``view3D-KINC-tidy.py``:  a Pyhon Dash application that uses Plotly to create 3D interactive viewer of the network.  It allows you to explore the co-expression relationships of each edge.  It requires a variety of python libraries to be installed prior to use.  This script is **only applicable to the GMM approach** of network.

Other Useful Functions
  - ``export-emx``: Export a gene expression matrix, or GEM, from the binary format used by KINC into a tab-delimited format compatible with other tools.
  - ``export-cmx``: Export a correlation matrix (or similarity matrix) from the binary format that KINC uses into a text-based version.
  - ``import-cmx``: Import a correlation matrix (or similarity matrix) that was previously exported by KINC.
