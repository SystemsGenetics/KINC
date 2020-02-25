Prepare Input Data
====================
Create the GEM
--------------
Before using KINC, you must have a valid Gene Expression Matrix (GEM) file. Please see the :doc:`data_overview` section for a description of the GEM file format. The GEM file can be created however is easiest for you.  If you are working with raw RNA-seq FASTQ files or would like to use RNA-seq data from the `NCBI sequence read archive <https://www.ncbi.nlm.nih.gov/sra>`_ (SRA) then you may consider using `GEMmaker <https://gemmaker.readthedocs.io/en/latest/>`_. GEMmaker is a sister tool of KINC, and is a Nextflow workflow designed to process large-scale RNA-seq datasets yielding a GEM file. It can create GEMs using tools such as `Hisat2 <https://ccb.jhu.edu/software/hisat2/index.shtml>`_, `Kallisto <https://pachterlab.github.io/kallisto/>`_ or `Salmon <https://combine-lab.github.io/salmon/>`_.

.. figure:: images/GEMmaker-logo-sm.png
   :alt: GEMmaker Logo

Normalizing the GEM
-------------------
If the GEM contains raw FPKM, RPKM or TPM values directly from the gene quantification tool, you must `log2` transform the data. This can be performed easily using R or Python. Sometimes, the samples must also undergo normalization.  Quantile normalization is often performed, but the most appropriate form of normalization is an open research question.

Create the Annotation Matrix
----------------------------
To generate condition-specific subgraphs you should also prepare a sample annotation matrix. Please see the :doc:`data_overview` section for a description of the file format. In short, this file should contain all of the experimental information about each sample as well as any phenotypic data about the individuals from which the samples were taken.
