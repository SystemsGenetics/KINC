Prepare the GEM
===============
Create the GEM
--------------
Before using KINC, you must have a valid Gene Expression Matrix (GEM) file.  The GEM file is a tab-delimited text file containing an `n` x `m` matrix of gene expression values where rows (or lines in the file) are genes, columns are samples and  elements are the expression values for a gene in a given sample.

The GEM may have a header line. But the header line must only contain the sample names. Each subsequent line will always start with the gene name, followed by the expression values.  Therefore, the header line, if present, should always have one less element then every other row.

.. note::

  The header line, if present, should always have one less element then every other row.

The GEM file can be created however is easiest for you.  If you are working with raw RNA-seq FASTQ files or would like to use RNA-seq data from the `NCBI sequence read archive <https://www.ncbi.nlm.nih.gov/sra>`_ (SRA) then you may consider using `GEMmaker <https://gemmaker.readthedocs.io/en/latest/>`_. GEMmaker is a sister tool of KINC, and is a Nextflow workflow designed to process large-scale RNA-seq datasets yielding a GEM file. It can create GEMs using tools such as `Hisat2 <https://ccb.jhu.edu/software/hisat2/index.shtml>`_, `Kallisto <https://pachterlab.github.io/kallisto/>`_ or `Salmon <https://combine-lab.github.io/salmon/>`_.

.. figure:: images/GEMmaker-logo-sm.png
   :alt: GEMmaker Logo

Normalizing the GEM
-------------------
If the GEM contains raw FPKM, RPKM or TPM values directly from the gene quantification tool, you must `log2` transform the data. This can be performed easily using R or Python. Sometimes, the samples must also undergo normalization.  Quantile normalization is often performed, but the most appropriate form of normalization is an open research question.
