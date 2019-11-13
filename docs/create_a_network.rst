How to Create a Network
=======================
This section describes the two approaches for network construction as well as some general considerations.  For the examples shown below we will assume the input GEM is derived from the SRA Project with ID `PRJNA301554 <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA301554/>`_. This dataset consists of 475 RNA-seq Illumina samples of rice grown under control, heat, drought, heat recovery and drought recover conditions.  It measures expression across four different genotypes of rice (from two subspecies) over a series of time points.  This dataset was selected because of its large size and multiple variables.

We will assume that these 475 samples have been processed using `GEMmaker <https://gemmaker.readthedocs.io/en/latest/>`_ and the resulting GEM has been `log2` transformed.

Before Getting Started
----------------------

Traditional or GMM Approach?
````````````````````````````
Before proceedng you should identify if a traditional or GMM network is most appropriate. The following considerations can help with that decision.

Traditional Approach
::::::::::::::::::::

Advantages
  - Executes very fast.
  - Biological signal can be found at higher correlation values and when sample bias is skewed towards the question of interest

Disadvantages
  - Includes false edges but also misses edges due to improper application of correlation tests.

GMM Approach
::::::::::::

Advantages
  - Limits correlation test bias by using GMMs to ensure test assumptions are met.
  - Can identify condition-specific pairwise expression.
  - Reduces the effect of sample bias when multiple conditions are present in the input GEM dataset.
  - Can execute very fast on a single machine for small GEMs.

Disadvantages
  - May require access to HPC with multiple GPUs for very large GEMs.
  - May not be useful if the sample metadata does not have qualitative (categorical) conditions.


.. _samples-needed-reference-label:

How Many Samples are Needed?
````````````````````````````

Networks can be created with very few samples if need be, but the power of the network will diminish greatly.  For Traditional networks, you can manually perform a power analysis before network construction to identify what correlation value (i.e. effect size) you must not go below in thresholding in order to limit false edges (assuming correlation assumptions are met, which they are not in the traditional approach). The ``pwr.r.test`` function of the statistical programming language R can do this, and there are `online calculators <http://www.sample-size.net/correlation-sample-size/>`_ as well.

For example, the minimum number of samples required to meet the criteria of for a significance value of 0.001, a power value of 0.8 and a minimum correlation threshold of 0.5 is 60 samples. If we raise the minimum threshold to 0.7 we need at least 21 samples.  Only 11 samples are needed for a threshold limit of 0.9.  If we only had 11 samples we should not allow a correlation threshold below 0.9.

If you are using the GMM approach and you wish to find condition-specific subgraphs for a qualitative condition, such as for genes underlying response to a treatment (e.g.heat, drought, control, etc.) you must ensure that you have sufficient samples for each category.  Suppose you only had 10 samples per treatment, you would expect to find clusters of approximately size 10 for each treatment, and this would require a minimum correlation threshold of 0.9. You can remove edges whose correlation dips below the limit using the ``corrpower`` function. You can set the minimum cluster size when executing the ``similarity`` function.

.. note::

  The number of samples will dictate the quantity and size of the final network.  With few samples sizes there is little chance of finding weakly correlated, but perhaps meaningful edges.

Do Replicates Matter?
`````````````````````
Unlike with DGE analysis, where multiple replicates are necessary to establish a statistical difference between two conditions, for GCNs a relationship between genes can be established using all data for that gene.  Therefore, It is not clear what contribution replicates make in the network.  Technical replicates are probably not useful.  Biological replicates may be useful but the benefit or bias from using replicates is not yet known.

Which Correlation Method to Choose?
```````````````````````````````````
KINC currently provides two correlation methods:  Pearson and Spearman.  Pearson is meant for linear relationships with no outliers that have equal variance.  Spearman is less affected by outliers and non-linearity so long as the relationship is monotonically increasing.  It is therefore recommended to use Spearman as not all relationships are linear.  However, Spearman correlation is less effective with small samples.  Empirically, Spearman tends to suffer when sample sizes dip below 20 to 30.  If you wish to identify condition-specific edges where each category (e.g. heat, drought, control) has fewer than 20 to 30 samples you should consider using Pearson.

Traditional Approach
--------------------
You can construct a traditional network on a stand-alone workstation using either ``kinc`` or ``qkinc``.  Using the 475-sample rice dataset described above, the following steps show how to create a traditional network using the command-line. The arguments shown in the command-line examples below correspond directly to fields in the KINC GUI.

Step 1: Import the GEM
``````````````````````
The first step is to import the GEM into a binary format suitable for KINC. The ``import-emx`` function of KINC does this:

.. code:: bash

  kinc run import-emx \
    --input "rice_heat_drought.GEM.FPKM.filtered.txt" \
    --output "rice_heat_drought.GEM.FPKM.filtered.emx" \
    --nan "NA" \
    --samples 0

In the example code above the GEM file is provided to the ``--input`` argument and the name of an output EMX file is provided using the ``--output`` argument.  In the example above, the ``--nan`` argument indicates that the file uses ``"NA"`` to represent missing values. This value should be set to whatever indicates missing values. This could be ``"0.0"``, ``"-Inf"``, etc. and the GEM file has a header describing each column so the number of samples provided to the ``--samples`` argument is set to 0. If the file did not have a header the number of samples would need to be provided.

Step 2: Perform Correlation Analysis
````````````````````````````````````
Construction of a similarity matrix (or correlation matrix) is the second step. Here KINC performs pairwise comparison of every gene with every other gene using either Spearman or Pearson correlation.  The ``similarity`` function of KINC does this:

.. code:: bash

  kinc run similarity \
    --input "rice_heat_drought.GEM.FPKM.filtered.emx" \
    --ccm "rice_heat_drought.GEM.FPKM.filtered.traditional.ccm" \
    --cmx "rice_heat_drought.GEM.FPKM.filtered.traditional.cmx" \
    --clusmethod "none" \
    --corrmethod "spearman" \
    --minsamp 30 \
    --minexpr -inf \
    --mincorr 0.5 \
    --maxcorr 1

Here the EMX file created in the first step is provided using the ``--emx`` argument and the names of two output files are provided using the ``--cmx`` and ``--ccm`` arguments. These are the correlation matrix and clustering matrix  respectively.  Because we are using the traditional approach, the ``--clusmethod`` argument is set to ``"none"``.  The correlation method is set to use Spearman, and the minimum number of samples required to perform correlation is set to 30 using the ``--minsamp`` argument. Any gene pairs where one gene has fewer that ``--minsamp`` samples will be excluded.  This will exclude genes that have missing values in samples that causes the number of samples to dip below this level.  The ``--minsamp`` argument should be set equal to or lower than the number of samples present in the origin GEM input file and higher than an expected level of missigness (e.g. 10% missing values allowed).  The ``--minexp`` argument isset to negative infinity (``-inf``) to indicate there is no limit on the minimum expression value.  If we wanted to exclude samples whose log2 expression values dipped below 0.2, for instance, we could do so with this argument.  To keep the output files relatively small, we will exclude all correlation values below 0.5 using the ``--mincorr`` argument.  Sometimes errors occur in data collection or quantification yielding high numbers of perfectly correlated genes!  We can limit that by excluding perfectly correlated genes by lowering the ``--maxcorr`` argument. In practice we leave this as 1 for the first time we create the network, if we fail to find a proper threshold in a later step then one cause may be large numbers of perfectly correlated genes.

Step 3: Thresholding
````````````````````
There are four ways KINC can determine a threhsold for a network: power-law, Random Matrix Theory (RMT), condition-specific and `ad hoc`.

.. _rmt-reference-label:

Method 1: Using RMT to Threshold
::::::::::::::::::::::::::::::::

The following command-line provides an example for RMT thresholding of the example 475-rice sample data:

.. note::

  RMT works best for traditional networks.

.. code:: bash

  kinc run rmt \
    --input "rice_heat_drought.GEM.FPKM.filtered.traditional.cmx" \
    --log "rice_heat_drought.GEM.FPKM.filtered.traditional.rmt.log" \
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

The above command provides the correlation matrix (CMX) using the ``--input`` arugment, and the name of a log file, using the ``--log`` argument  where the results of chi-square testing is stored.  The RMT method will successively walk through all correlation values, in decreasing order from ``--tstart`` to ``--tstop``, using a step of ``--tstep``, and builds a new similarity matrix to test if the Nearest Neighbor Spacing Distribution (NNSD) of the Eigenvalues of that matrix appears Poisson.  A spline curve is fit to the NNSD if the ``--spline`` argument is ``TRUE`` (recommended) and random points along the line are selected to determine if the distribution appears Poisson.  This random selection will occur repeatedly by selecting a random set of ``--minpace`` numbers and increasing that on successive iterations to ``--maxpace``.  A Chi-square test is performed for each of these random selections and the result is averaged for each correlation value.  The ``--bins`` is the number of bins in the NNSD histogram and `1 - bins` indicates how many degrees of freedom the Chi-square test will have. In practice, a Chi-square value of 100 indicates that the correlation value begins to not look Poisson. The RMT approach will continue after seeing a Chi-square value of 100 until it sees one at the 200 at which point it stops.  It seeks past 100 to ensure it does not get trapped in a local minimum.

.. note::

  It is best to leave all options as default unless you know how to tweak the RMT process.

Once completed, you can determine the best threshold for the network by opening the logfile specified by the ``--log`` argument, and looking at the end of the file.  The threshold is listed on the last line of the file and should be used for extracing the network in step 4.

If the input GEM is especially noisy, the RMT method will fail to find a threshold. As it continues to search through decreasing correlation values, the time required to generate the eigenvalues dramatically increases and it may appear that RMT never completes.  To determine if this is the case, examine the log file. If you see the average correlation beyond 200 then this has occurred.  See the :doc:`troubleshooting` section to explore alternative methods.


Method 2: Using the Power-law Threshold
:::::::::::::::::::::::::::::::::::::::
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

  While the power-law threshold is useful to help identify scale-free behavior, it does not that the network is modular and hierarchical.

Method 3: Applying a Condition-Specific Filter
::::::::::::::::::::::::::::::::::::::::::::::
The condition-specific thresholding approach uses an annotation matrix that contains metadata about the samples such as the experimental conditions or phenotypes.  The approach to perform condition-specific thresholding is the same as for the GMM approach. Please refer to the :ref:`csfilter-reference-label` section for details about using condition-specific filters for either traditional or GMM networks.

.. warning::

  Condition-specific thresholding only works with traditional networks when the metadata in the annotation matrix is quantitative.

Method 2: Using an `Ad Hoc` Approach
::::::::::::::::::::::::::::::::::::
An `ad hoc` threshold does not use an anlytical approach to determine a threshold. Instead, the researcher selects a reasonable threshold. For example, this could be the minimum correlation that selects the top 1000 relationships, or yields a network that has desired size or communities.  These types of thresholding approaches have been used for peer-reviewed published networks but users should be cautious when using this approach.

Step 4: Extracting the Network File
```````````````````````````````````
How ever you have chosen to threshold the network, either with RMT or Power-law, or some `ad-hoc` approach, you will have a minimum correlation value.  This value can be used to extract any pairwise comparison between genes in the correlation matrix file (CMX) that are above the absolute value of the minimum correlation. These become edges in the final network.  The ``extract`` function of KINC will do this:

.. code:: bash

  kinc run extract \
     --emx "rice_heat_drought.GEM.FPKM.filtered.emx" \
     --ccm "rice_heat_drought.GEM.FPKM.filtered.traditional.ccm" \
     --cmx "rice_heat_drought.GEM.FPKM.filtered.traditional.cmx" \
     --format "text" \
     --output "rice_heat_drought.GEM.FPKM.filtered.traditional.gcn.txt" \
     --mincorr 0.892001 \
     --maxcorr 1

As in previous steps, the ``--emx``, ``--cmx`` and ``--ccm`` arguments provide the exrpession matrix, correlation and clustering matricies. The threshold is provided to the ``--mincorr`` argument.  Additinally, if you would like to exclude high correlations (such as perfect correlations), you can do so with the ``--maxcorr`` argument. You should only need to change the ``--maxcorr`` argument if it was determined that there is error in the data resulting in an inordinate number of high correlations.  The ``--format`` argument can be ``text``, ``minimal`` or ``graphml``. The ``text`` format currently contains the most data. It is easily imported into Cytoscape or R for other analyses and visualizations. The ``minimal`` format simply contains the list of edges with only the two genes and the correlation value. The ``graphml`` format provides the same information as the ``minimal`` format but using the `GraphML <http://graphml.graphdrawing.org/>`_ file format.

See the :ref:`plain-text-reference-label`  section for specific details about these files.

GMM approach
------------
Here we perform network construction using the Gaussian Mixture Model (GMM) appraoch.  With this approach, each pair-wise comparision of every two genes undergoes a cluster identification analysis using GMMs. This approach can identify clusters, or groups, of samples that have similar but distinct ranges of expression. The underlying hypothesis is that when clusters appear, they represent condition-specific gene expression.  Clusters that are identified in gene pairs are correlated independently and each cluster has the potential to become a separate edge in the network.  Because we know the samples that are present in each cluster, KINC uses a hypergeometric test to compare categorical data about samples with cluster membership, and regression analysis to compare qualitative and ordinal data. Condition-specific thresholding can be performed on the `p`-values and `r`-squared values of those test to generate condition-specific subgraphs.

.. note::

  The GMMs approach requires a tab-delimited annotation matrix file (AMX) that contains metadata about samples where columns are feature that contain experimental condition information or phenotype data.

Step 1: Import the GEM
``````````````````````
.. code:: bash

  kinc run import-emx \
    --input "rice_heat_drought.GEM.FPKM.filtered.txt" \
    --output "rice_heat_drought.GEM.FPKM.filtered.emx" \
    --nan "NA" \
    --samples 0

In the code above the GEM file is provided to the ``import-emx`` function and the name of an output EMX file is provided.  The file uses "NA" to indicate missing values and  it has a header so the number of samples is set to .

Step 2: Perform GMM + Correlation Analysis
``````````````````````````````````````````
The second step is to use GMM to identify clusters and then perform correlation analysis on each cluster.

.. code:: bash

  kinc run similarity \
    --input "rice_heat_drought.GEM.FPKM.filtered.emx" \
    --ccm "rice_heat_drought.GEM.FPKM.filtered.ccm" \
    --cmx "rice_heat_drought.GEM.FPKM.filtered.cmx" \
    --clusmethod "gmm" \
    --corrmethod "spearman" \
    --minexpr -inf \
    --minsamp 25 \
    --minclus 1 \
    --maxclus 5 \
    --crit "ICL" \
    --preout TRUE \
    --postout TRUE \
    --mincorr 0.5 \
    --maxcorr 1

Here the EMX file created in the first step is provided, and the names of the two output (CCM and CMX) files are provided.  Because we are using the GMM approach, the ``--clusmethod`` argument is set to ``"gmm"``.  The correlation method is set to use Spearman.  Other argument specific to the GMM appraoch include ``--crit``, ``--maxclus``, ``-minclus``, ``--preout``, and ``--postout``. These have the following meaning:

-  ``--crit``: This is the criterion to select a clustering model. This should remain as ``ICL`` unless a higher number of modules per pair is desired and can be set to 'BIC'.
- ``--minclus``: The minimum number of clusters that can be found per gene pair.  Unless you are specifically looking for genes with multi-modal expression this should remain s 1.
- ``--maxclus``: The maximum number of clusters that can be found per gene pair.
- ``--preout``: Set to TRUE to turn on removal of outliers prior to GMM clustering. FALSE otherwise.
- ``--postout``:  Set to TRUE to remove outliers that may be present in GMM clusters. FALSE  otherwise.


The ``--minexp`` argument isset to negative infinity (``-inf``) to indicate there is no limit on the minimum expression value.  If we wanted to exclude samples whose log2 expression values dipped below 0.2, for instance, we could do so.  To keep the output files relatively small, we will exclude all correlation values below 0.5 using the ``--mincorr`` argument.

Sometimes errors occur in data collection or quantification yielding high numbers of perfectly correlated genes!  We can limit that by excluding perfectly correlated genes by lowering the ``--maxcorr`` argument. In practice we leave this as 1 for the first time we create the network.


Step 3: Filter Low-Powered Edges
````````````````````````````````
As discussed in the :ref:`samples-needed-reference-label` section above, the power of a correlation analysis is dependent on the number of samples in the test.  Unlike the traditional approach, where a power analysis can indicate the minimum correlation threshold below which you should not drop, a power-analysis for the GMM approach must be applied to each cluster separately.  The ``corrpower`` function does this and removes underpowered clusters from the matricies. For example:

.. code:: bash

  kinc run corrpower \
  --ccm-in "rice_heat_drought.GEM.FPKM.filtered.ccm" \
  --cmx-in "rice_heat_drought.GEM.FPKM.filtered.cmx" \
  --ccm-out "rice_heat_drought.GEM.FPKM.filtered.paf.ccm" \
  --cmx-out "rice_heat_drought.GEM.FPKM.filtered.paf.cmx" \
  --alpha 0.001 \
  --power 0.8

As shown above, the power and signficance criteria are set with the ``--power`` and ``--alpha`` arguments respectively.  An ``alpha`` setting of ``0.001`` indicates that we want to limit the Type I error (false positives) to a signicance level of 0.001.  The Power uses the formula 1-`Beta` where `Beta` is the probability of a Type II error (false negative) occuring.  A ``--power`` setting of 0.8 indicates that we are comfortable with a 20% false negative rate. There is no rule for how to set these.  Set them to the levels of noise you are comfortable with.

.. note::

  Remember, to find edges in the nework associated with categorical features, you must have enough samples with the given category in order to find a cluster an then to have sufficent power. The ``--minsamp `` argument in the ``similarity`` step sets the smallest allowable cluster size.

Step 4: Thresholding
````````````````````
For the GMM aproach there are several options for thresholding: Random Matrix Theory (RMT) or condition-specific.

Method 1: Using RMT to Threshold
::::::::::::::::::::::::::::::::

You can use RMT for identifying a threshold for the GMM approach.  For this you should adjust the ``--reduction`` argument to specify one of: ``first``, ``mincorr``, ``maxcorr`` or ``random``.  This will select the cluster that  appears first, has the minimum correlation value, maximum correlation value or a random cluster, respectively.  However, RMT cannot be used for identifying condition-specific subgraphs and results in a traditional style network even if the GMM approach was used. To use RMT please see the :ref:`rmt-reference-label` section.


.. _csfilter-reference-label:

Method 2: Applying a Condition-Specific Filter
::::::::::::::::::::::::::::::::::::::::::::::
Condition-specific filtering is performed using the ``cond-test`` function of KINC. It requires an annotation matrix containing metadata about the RNA-seq samples. It performs a hypergeometric test for categorical features and regression analysis for quantitative features and assignes `p`-values and `r`-squared values, as appropriate, to each edge in the network. The following shows an example:

.. code:: bash

  kinc run cond-test \
    --emx "rice_heat_drought.GEM.FPKM.filtered.emx" \
    --ccm "rice_heat_drought.GEM.FPKM.filtered.paf.ccm" \
    --cmx "rice_heat_drought.GEM.FPKM.filtered.paf.cmx" \
    --amx "../../01-input_data/rice_heat_drought/PRJNA301554.hydroponic.sample_annotations.filtered.txt"   \
    --output "rice_heat_drought.GEM.FPKM.filtered.paf.csm" \
    --feat-tests "Subspecies,Treatment,GTAbbr" \
    --feat-types "Subspecies:categorical,Treatment:categorical:GTAbbr:categorical"

Here, the ``--emx``, ``--ccm``, and ``--cmx`` arguments provide the usual expression matrix, cluster matrix and correlation matrix respectively.  The ``--amx`` argument specifies the :ref:`amx-reference-label`.  The name of new condition-specific matrix, that will contain the results of the tests is set using the  ``--output`` argument.

Finally, it may not be desired to test all of the metadata features (i.e. columns) from the annotation matrix.  Using the ``feat-tests`` argument you can specify a comma-separated list (without spaces) of the names of the columns in the annotation matrix file that should be tested.  These can be either categorical, quantitative or ordinal.  KINC will do its best to determine the top of data in each column, but you can override the type using the ``--feat-types`` argument and specifying the type by separating with a colon.

Unlike with other thresholding methods, you do not get a minimal correlation value. Instead you can set limits on the `p`-values and `r`-squared values in the network extraction step.

Step 5: Extract Condition-Specific Subgraphs
````````````````````````````````````````````
When using the GMM approach, the goal is to identiy condition-specific subgraphs. These are subsets of a larger "unseen" network that are specific to a given condition.  As with the traditional approach, the ``extract`` function of KINC will do this:

.. code:: bash

  kinc run extract \
    --emx "rice_heat_drought.GEM.FPKM.filtered.emx" \
    --ccm "rice_heat_drought.GEM.FPKM.filtered.paf.ccm" \
    --cmx "rice_heat_drought.GEM.FPKM.filtered.paf.cmx" \
    --csm "rice_heat_drought.GEM.FPKM.filtered.paf.csm" \
    --format "text" \
    --output "rice_heat_drought.GEM.FPKM.filtered.th0.5.cs1e-3.gcn.txt" \
    --mincorr 0.80 \
    --maxcorr 1 \
    --filter-pvalue "1e-3"
    --filter-rsquare "0.3"



As in previous steps, the ``--emx``, ``--cmx``, ``--ccm`` and ``--csm`` arguments provide the exrpession matrix, correlation,  clustering matrix and the new condition-specific matrix. A threshold is provided to the ``--mincorr`` argument typically as a lower-bound. No edges with absolute correlation values below this value will be extracted.   Additinally, if you would like to exclude high correlations (such as perfect correlations), you can do so with the ``--maxcorr`` argument. You should only need to change the ``--maxcorr`` argument if it was determined that there is error in the data resulting in an inordinate number of high correlations.  To limit the size of the condition-specific subgraphs you should then set the ``--filter-pvalue`` and ``--filter-rsquare`` values to lower-bounds for signficant p-values and meaningful r-square values from test.  The r-square values are only present for quantitative features where the regression test was performed.  The p-value in this case indicates how well the data follows a trend and the r-square indicates how much of the variation the trend line accounts for.  Ideally, low p-values and high r-squre are desired. However, there are no rules for the best setting, but choose settings that provide a signficance level you are comfortable with.

Finally, the ``--format`` argument can be ``text``, ``minimal`` or ``graphml``. The ``text`` format currently contains the most data. It is easily imported into Cytoscape or R for other analyses and visualizations. The ``minimal`` format simply contains the list of edges with only the two genes and the correlation value. The ``graphml`` format provides the same information as the ``minimal`` format but using the `GraphML <http://graphml.graphdrawing.org/>`_ file format.

See the :ref:`plain-text-reference-label`  section for specific details about these files.


Step 6: Remove Edges Due to Collinearity
````````````````````````````````````````
This last step must be performed in R using the KINC.R package.  KINC is an actively developed software project and functions are often implemented in R before beinc moved to the faster C++ based KINC software.  Currently, the function for removing edges due to collinearity is in the KINC.R supplemental package.  To use this package, you must first have all dependencies installed these include the following CRAN and Bioconductor packages:

CRAN

- ggplot2
- igraph
- linkcomm
- Rmixmod
- pwr
- dplyr
- reshape2
- grid

Bioconductor

- qvalue

To install KINC.R, follow the instructions on the `KINC.R Github repository <https://github.com/SystemsGenetics/KINC.R>`_.

Once KINC.R and all dependencies you can execute the following:

Frist, import the exrpession and annotation matricies (the same used for KINC):

.. code:: r

  ematrix = read.table('rice_heat_drought.GEM.FPKM.filtered.txt', header=TRUE, sep="\t")
  osa = loadSampleAnnotations('PRJNA301554.hydroponic.sample_annotations.filtered.txt')
  net = loadKINCNetwork('rice_heat_drought.GEM.FPKM.filtered-th0.85-p1e-10.r0.3-gcn.txt')

Next we should filter edges that may be biased due to collinearity:

.. code:: r

  net2 = filterBiasedEdges(net, ematrix, th=1e-3)

Now that the filter has been applied we can save the final network file:
saveKINCNetwork(net2, 'rice_heat_drought.GEM.FPKM.filtered-th0.85-p1e-10.r0.3-b1e-3.gcn.txt')

..warning::

  It is possible that condition-specific networks that are brought into KINC.R can be extremely large. Such is the case when time series are present in the data resulting in many collinear relationships due to circadian changes.  To avoid overrunning memory in R, you may want to filter your next during the ``extract`` phase using a higher p-value threshold or minimum correlation value to limit the size of the network.
