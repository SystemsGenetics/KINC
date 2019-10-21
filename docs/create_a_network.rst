Create a Network
================
This section describes the two approaches for network construction as well as some general considerations.  For the examples shown below we will assume the input GEM is derived from the SRA Project with ID `PRJNA301554 <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA301554/>`_. This dataset consists of 475 RNA-seq Illumina samples of rice grown under control, heat, drought, heat recovery and drought recover conditions.  It measures expression across four different genotypes of rice (from two subspecies) over a series of time points.  This dataset was selected because of its large size and multiple variables.

We will assume that these 475 samples have been processed using `GEMmaker <https://gemmaker.readthedocs.io/en/latest/>`_ and the resulting GEM has been `log2` transformed.

How Many Samples are Needed?
----------------------------

Which Correlation Method to Choose?
-----------------------------------



Traditional or GMM Approach?
----------------------------
he traditional approach to network construction performs a single pairwise correlation test between all genes.  The advantage is that it executes relatively quickly. The disadvantage is that it does not ensure that the assumptions of the correlation test are met, leading to false edges and missing edges. Traditional networks have been useful for finding some biological signal for very strongly correlated relationships and when sample bias is skewed towards the biological question of interest.



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
