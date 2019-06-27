Analytics
=========

KINC consists of several sub-commands called "analytics", each of which performs a single step in the network construction process. Each analytic is described here.

To view the help text for an analytic, run ``kinc help run <analytic>``. The help text lists all of the command-line options for the analytic.

Import Expression Matrix
------------------------

The Import Expression Matrix analytic takes a plain-text expression matrix and converts it to an EMX object. Elements which have the given NAN token ("NA" by default) are read in as NAN. If the sample names are not in the input file, the user must provide the number of samples to the analytic, and the samples will be given integer names.

.. code:: bash

   kinc run import-emx --input GEM.emx.txt --output GEM.emx

Export Expression Matrix
------------------------

The Export Expression Matrix analytic takes an EMX object and converts it to a plain-text expression matrix. Elements which are NAN in the expression matrix are written as the given NAN token ("NA" by default).

.. code:: bash

   kinc run export-emx --input GEM.emx --output GEM.emx.txt

Similarity
----------

The Similarity analytic takes an EMX object and computes a similarity matrix, where each element is a similarity measure of two genes in the expression matrix. The similarity of each pair is computed using a correlation measure. The similarity matrix can also have multiple modes within a pair; these modes can be optionally computed using a clustering method. The similarity matrix produced by this analytic consists of two data objects: a CMX object containing the pairwise correlations, and a CCM object containing sample masks of the pairwise clusters. If clustering is not used, an empty CCM object is created. This analytic can also perform pairwise outlier removal before and after clustering, if clustering is used.

This analytic can use MPI and it has both CPU and GPU implementations, as the pairwise clustering significantly increases the amount of computations required for a large expression matrix.

.. code:: bash

   kinc run similarity --input GEM.emx --ccm GEM.ccm --cmx GEM.cmx

Export Correlation Matrix
-------------------------

The Export Correlation Matrix analytic takes two data objects, a CMX object and a CCM object, and writes a plain-text correlation matrix.

.. code:: bash

   kinc run export-cmx --emx GEM.emx --ccm GEM.ccm --cmx GEM.cmx --output GEM.cmx.txt

Import Correlation Matrix
-------------------------

The Import Correlation Matrix analytic takes a plain-text correlation matrix and produces two data objects: a CMX object containing the pairwise correlations, and a CCM object containing the sample masks for each pairwise cluster. There are several fields which are not represented in the input file and therefore must be specified manually, including the gene size, sample size, max cluster size, and correlation name.

.. code:: bash

   kinc run import-cmx --input GEM.cmx.txt --ccm GEM.ccm --cmx GEM.cmx

Correlation Power Filtering
---------------------------

The Correlation Power Filtering analytic takes a CMX object and CCM object and filters out clusters that do not pass specific tests. This analytic can use MPI.

.. code:: bash

   kinc run corrpower \
      --ccm-in GEM.ccm \
      --cmx-in GEM.cmx \
      --ccm-out GEM.corrpower.ccm \
      --cmx-out GEM.corrpower.cmx

Power-law Thresholding
----------------------

The Power-law Thresholding analytic takes a CMX object and attempts to find a threshold which, when applied to the correlation matrix, produces a scale-free network. Each thresholded network is evaluated by comparing the degree distribution of the network to a power-law distribution. This process is repeated at each threshold step from the starting threshold to the stopping threshold.

.. code:: bash

   kinc run powerlaw --input GEM.cmx --log GEM.powerlaw.txt

RMT Thresholding
----------------

The RMT Thresholding analytic takes a CMX object and attempts to find a threshold which, when applied to the correlation matrix, produces a non-random network. This analytic uses Random Matrix Theory (RMT), which involves computing the eigenvalues of a thresholded correlation matrix, computing the nearest-neighbor spacing distribution (NNSD) of the eigenvalues, and comparing the distribution to a Poisson distribution using a chi-squared test. This process is repeated at each threshold step from the starting threshold; as the threshold decreases, the NNSD transitions from a Poisson distribution to a Gaussian orthogonal ensemble (GOE) distribution, which causes the chi-squared value to increase sharply. The final threshold is chosen as the lowest threshold whose chi-squared value was below the critical value.

.. code:: bash

   kinc run rmt --input GEM.cmx --log GEM.rmt.txt

Extract
-------

The Extract analytic takes two data objects, a CMX object and a CCM object, and extracts a co-expression network by applying a correlation threshold. The network file can use the full plain-text format, the "minimal" plain-text format, or the GraphML format.

.. code:: bash

   kinc run extract --emx GEM.emx --ccm GEM.ccm --cmx GEM.cmx --output GEM.coexpnet.txt
