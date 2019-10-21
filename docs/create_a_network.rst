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
Networks can be created with very few samples if need be, but the power of the network will diminish greatly.  For Traditional networks, you can manually perform a power analysis before network construction to identify what correlation value (i.e. effect size) you must not go below in thresholding in order to limit false edges (assuming correlation assumptions are met, which they are not in the traditional approach). The ``pwr.r.test`` function of the statistical programming language R can do this, and there are `online calculators <http://www.sample-size.net/correlation-sample-size/>`_ as well.

For example, the minimum number of samples for an alpha (Type I error) value of 0.001, a beta (Power: 1 - Type II error) value of 0.8 and a minimum correlation threshold of 0.5 is 60 samples. If we raise the minimum threshold to 0.7 we need at least 21 samples.  Only 11 samples are needed for a threshold limit of 0.9.  If we only had 11 samples we should not allow a correlation threshold below 0.9.

If you are using the GMM approach and you wish to find condition-specific subgraphs for a qualitative condition, such as for genes underlying response to a treatment (e.g.heat, drought, control, etc.) you must ensure that you have sufficient samples for each category.  Suppose you only had 10 samples per treatment, you would expect to find clusters of approximately size 10 for each treatment, and this would require a minimum correlation threshold of 0.9. You can remove edges whose correlation dips below the limit using the ``corrpower`` function. You can set the minimum cluster size when executing the ``similarity`` function.

.. note::

  The number of samples will dictate the quantity and size of the final network.  With few samples sizes there is little chance of finding weakly correlated, but perhaps meaningful edges.

Do Replicates Matter?
`````````````````````
Unlike with DGE analysis, where multiple replicates are necessary to establish a statistical difference between two conditions, for GCNs a relationship between genes can be established using all data for that gene.  Therefore, It is not clear what contribution replicates make in the network.  Technical replicates are probably not useful.  Biological replicates may be useful but the benefit or bias from using replicates is not yet known.

Which Correlation Method to Choose?
```````````````````````````````````
KINC currently provides two correlation methods:  Pearson and Spearman.  Pearson is meant for linear relationships with no outliers that have equal variance.  Spearman is less affected by outliers and non-linearity so long as the relationship is monotonically increasing.  It is therefore recommended to use Spearman as not all relationships are linear.  However, Spearman correlation is less effective with small samples.  Empirically, Spearman tends to suffer when sample sizes dip below 20 to 30.  If you wish to identify condition-specific edges where each category (e.g. heat, drought,control) has fewer than 20 to 30 samples you should consider using Pearson.

Traditional Approach
--------------------
You can construct a traditional network on a stand-alone workstation using either ``kinc`` or ``qkinc``.  Using the 475-sample rice dataset described above, the following steps show how to create a traditional network using the command-line. The arguments shown in the command-line examples below correspond directly to fields in the KINC GUI.

Step 1: Import the GEM
``````````````````````
.. code:: bash

  kinc run import-emx \
    --input "rice_heat_drought.GEM.FPKM.filtered.txt" \
    --output "rice_heat_drought.GEM.FPKM.filtered.emx" \
    --nan "NA" \
    --samples 0

In the code above the GEM file is provided to the ``import-emx`` function and the name of an output EMX file is provided.  The file uses "NA" to indicate missing values and  it has a header so the number of samples is set to .

Step 2: Perform Correlation Analysis
````````````````````````````````````
.. code:: bash

  kinc run similarity \
    --input "rice_heat_drought.GEM.FPKM.filtered.emx" \
    --ccm "rice_heat_drought.GEM.FPKM.filtered.traditional.ccm" \
    --cmx "rice_heat_drought.GEM.FPKM.filtered.traditional.cmx" \
    --clusmethod "none" \
    --corrmethod "spearman" \
    --minexpr -inf \
    --mincorr 0.5 \
    --maxcorr 1

Here the EMX file created in the first step is provided, and the names of the two output (CCM and CMX) files are provided.  Because we are using the traditional approach, the ``--clusmethod`` argument is set to ``"none"``.  The correlation method is set to use Spearman and the ``--minexp`` argument isset to negative infinity (``-inf``) to indicate there is no limit on the minimum expression value.  If we wanted to exclude samples whose log2 expression values dipped below 0.2, for instance, we could do so.  To keep the output files relatively small, we will exclude all correlation values below 0.5 using the ``--mincorr`` argument.  Sometimes errors occur in data collection or quantification yielding high numbers of perfectly correlated genes!  We can limit that by excluding perfectly correlated genes by lowering the ``--maxcorr`` argument. In practice we leave this as 1 for the first time we create the network.

Step 3: Thresholding
````````````````````
There are three ways KINC can determine the threhsold for a network: Power-law analysis, Random Matrix Theory and Condition-Specific p-value limits.

Using RMT to Threshold
::::::::::::::::::::::
The following command-line executes RMT for the example 475-rice sample data:

.. note::

  RMT works best for traditional networks.

.. code:: bash

  kinc run rmt \
    --input "rice_heat_drought.GEM.FPKM.filtered.traditional.cmx" \
    --log "rice_heat_drought.GEM.FPKM.filtered.traditional.rmt.log" \
    --reduction "first" \
    --tstart 0.99 \
    --tstep 0.001 \
    --tstop 0.5 \
    --threads 1 \
    --epsilon 1e-6 \
    --mineigens 50 \
    --spline TRUE \
    --minpace 10 \
    --maxpace 40 \
    --bins 60

The above command provides the correlation matrix (CMX) file and the name of a log file where the results of the chi-square test are provided.  The RMT method will successively walk through all correlation values, in decreasing order from ``--tstart`` to ``--tstop``, using a step of ``--tstep``, and build a new similarity matrix to test if the Nearest Neighbor Spacing Distribution (NNSD) of the Eigenvalues appears Poisson.  A spline curve is fit to the NNSD if the ``--spline`` argument is ``TRUE`` (recommended) and random points along the line are selected to determine if the distribution appears Poisson.  This random selection will occur repeatedly by selecting a random set of ``--minpace`` numbers and increasing that on successive iterations to ``--maxpace``.  A Chi-square test is performed for each of these random selections and the result is averaged for each correlation value.  The ``--bins`` indicates how many degrees of freedom the Chi-square test will have. In practice, a Chi-square value of 100 indicates that the correlation value begins to not look Poisson. The RMT approach will continue after seeing a Chi-square value of 100 until it sees one at the 200 at which point it stops.  It seeks past 100 to ensure it does not get trapped in a local minimum.

In short, you can determine the best threshold for the network by opening the logfile specified by the ``--log`` argument, and looking at the end of the file.  The threshold that should be used for extracing the network is provided as the last number in the file.

In the input GEM is especially noisy, the RMT method will fail to find a threshold. As it continues to search through decreasing correlation values, the time required to generate eigenvalues dramatically increases and it may appear that RMT never completes.  To determine if this is the case examine the log file. If you see the average correlation beyond 200 then this has occurred.  See the :doc:`troubleshooting` section to explore alternative methods.

.. warning::

  You can use RMT for identifying a threshold for the GMM approach.  For this you should adjust the ``--reduction`` argument to specify one of: ``first``, ``mincorr``, ``maxcorr`` or ``random``.  This will select the cluster that  appears first, has the minimum correlation value, maximum correlation value or a random cluster, respectively.  However, RMT cannot be used for identifying condition-specific subgraphs and results in a traditional style network even if the GMM approach was used.

Using the Power-law Threshold
:::::::::::::::::::::::::::::
The Power-law function tests to see if the network, at successively decreasing correlation values follows a power-law which is a properly of scale-free network.  The power-law threshold can be used as an alternative to RMT when it fails to find a solution. The following example uses the power-law threshold for the example 475-rice sample data:

.. code:: bash

  kinc run powerlaw \
    --input "rice_heat_drought.GEM.FPKM.filtered.traditional.cmx" \
    --log "rice_heat_drought.GEM.FPKM.filtered.traditional.powerlaw.log" \
    --tstart 0.99 \
    --tstep 0.01 \
    --tstop 0.5

Here the correlation matrix (CMX) file is provided as well as a log file where details about the analysis are stored. The ``--tstart`` argument sets the starting correlation value and power-law calculations continue until the ``--tstop`` value is reached.

If function fails to find an threshold then see the :doc:`troubleshooting` section to explore alternative methods.

.. warning::

  While the power-law threshold is useful to help identify scale-free behavior, it cannot ensure that the network is modular and hierarchical.

Step 4: Extracting the Network File
```````````````````````````````````

GMM approach
------------
