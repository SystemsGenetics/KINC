Create a Network
================
This section describes the two approaches for network construction as well as some general considerations.  For the examples shown below we will assume the input GEM is derived from the SRA Project with ID `PRJNA301554 <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA301554/>`_. This dataset consists of 475 RNA-seq Illumina samples of rice grown under control, heat, drought, heat recovery and drought recover conditions.  It measures expression across four different genotypes of rice (from two subspecies) over a series of time points.  This dataset was selected because of its large size and multiple variables.

We will assume that these 475 samples have been processed using `GEMmaker <https://gemmaker.readthedocs.io/en/latest/>`_ and the resulting GEM has been `log2` transformed.

Before Getting Started
----------------------

Traditional or GMM Approach?
````````````````````````````
**Traditional Approach**

Advantages
  - Executes very fast.
  - Biological signal can be found at higher correlation values and when sample bias is skewed towards the question of interest

Disadvantages
  - Includes false edges but also misses edges due to improper application of correlation tests.

**GMM Approach**

Advantages
  - Limits correlation test bias by using GMMs to ensure test assumptions are met.
  - Can identify condition-specific pairwise expression.
  - Reduces the effect of sample bias when multiple conditions are present in the input GEM dataset.
  - Can execute very fast on a single machine for small GEMs.

Disadvantages
  - May require access to HPC with multiple GPUs for very large GEMs.
  - May not be useful if the experiment that generated the data does not have qualitative conditions. For example, an experiment to explore genes underlying differences in plant height across a range of genotypes would not benefit from GMMs as the phenotype is quantitative. A traditional approach would be better.

**KINC's Advantage for Both**

In both cases, when condition-specific thresholding is used, edges that are normally thrown out due to high correlation thresholds are not lost.  But edges that are correlated with qualitative conditions (via bionomial tests) of with quantitative conditions (via regression tests) are identified.

How Many Samples are Needed?
````````````````````````````
Networks can be created with very few samples if need be, but the power of the network will diminish greatly.  For Traditional networks, you can manually perform a power analysis before network construction to identify, what correlation value (i.e. effect size) you must not go below in thresholding in order to limit false edges (assuming correlation assumptions are met, which they are not in the traditional approach). The ``pwr.r.test`` function of the statistical programming language R can do this, and there are `online calculators <http://www.sample-size.net/correlation-sample-size/>`_ as well.

For example, the minimum number of samples for an alpha (Type I error) value of 0.001, a beta (Power: 1 - Type II error) value of 0.8 and a minimum correlation threshold of 0.5 is 60 samples. If we raise the minimum threshold to 0.7 we only need 21 samples.  Only 11 samples are needed for a threshold limit of 0.9.  If we only had 11 samples we should not allow a correlation threshold below 0.9.

If you are using the GMM approach and you wish to find condition-specific subgraphs for a qualitative condition, such as for genes underlying response to a treatment (e.g.heat, drought, control, etc.) you must ensure that you have sufficient samples for each category.  Suppose you only had 10 samples per treatment, you would expect to find clusters of approximately size 10 for each treatment, and this would require a minimum correlation threshold of 0.9. You can remove edges whose correlation dips below the limit using the ``corrpower`` function. You can set the minimum cluster size when executing the ``similarity`` function.

.. note::

  The number of samples will dictate the quantity and size of the final network.  With few samples sizes there is little chance of finding weakly correlated, but perhaps meaningful edges.

Do Replicates Matter?
`````````````````````
Unlike with DGE analysis, where multiple replicates are necessary to establish a statistical difference between two conditions, for GCNs a relationship between genes can be established using all data for that gene.  Therefore, It is not clear what contribution replicates make in the network.  Technical replicates are probably not useful.  Biological replicates may be useful but the benefit or bias from using replicates is not yet known.

Which Correlation Method to Choose?
```````````````````````````````````
KINC currently provides two correlation methods:  Pearson and Spearman.  Pearson is meant for linear relationships with no outliers that have equal variance.  Spearman is less affeccted by outliers and non-linearity so long as the relationship is monotonically increasing.  It is therefore recommended to use Spearman as not all relationships are linear.  However, Spearman correlation is less effective with small samples.  Empirically, Spearman tends to suffer when sample sizes dip below 30.

Traditional Approach
--------------------


GMM approach
------------

How to Choose a Minimum Cluster Size
````````````````````````````````````


To execute the functions described in the :doc:`functional_overview` section you will always

The easiest way to learn how to use KINC is to study the ``kinc.sh`` script, which can run the entire KINC workflow. Additionally, you can use the ``make-input-data.py`` script to generate a "fake" GEM with which to test KINC quickly:

.. code:: bash

   # generate fake GEM
   python scripts/make-input-data.py

   # run KINC
   scripts/kinc.sh serial 1 GEM.txt
