Quick Start
===========

This guide will get you started with KINC-Python in 5 minutes.

Installation
------------

Install KINC-Python using pip:

.. code-block:: bash

   pip install kinc-python

Basic Workflow
--------------

KINC-Python workflows typically follow these steps:

1. **Import** gene expression matrix from text file
2. **Compute similarity** (correlations) between gene pairs
3. **Extract** network edges to a file

Example with Sample Data
-------------------------

Let's walk through a complete example using synthetic data.

Step 1: Create Sample Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a simple gene expression matrix (GEM) file:

.. code-block:: bash

   cat > sample_data.txt << EOF
   	Sample1	Sample2	Sample3	Sample4	Sample5
   Gene1	1.2	2.3	1.5	2.1	1.8
   Gene2	3.4	4.1	3.8	4.2	3.9
   Gene3	0.5	0.7	0.6	0.8	0.7
   Gene4	2.1	1.9	2.0	1.8	2.2
   Gene5	4.5	4.3	4.7	4.1	4.6
   EOF

Step 2: Import to HDF5
~~~~~~~~~~~~~~~~~~~~~~

Convert the text file to HDF5 format:

.. code-block:: bash

   kinc import-emx --input sample_data.txt --output sample.emx.h5

Output:

.. code-block:: text

   Importing expression matrix from sample_data.txt...
     Genes: 5
     Samples: 5
   Expression matrix saved to sample.emx.h5

Step 3: Compute Similarity
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compute pairwise correlations using Pearson method:

.. code-block:: bash

   kinc similarity \
     --input sample.emx.h5 \
     --cmx sample.cmx.h5 \
     --ccm sample.ccm.h5 \
     --corrmethod pearson \
     --clusmethod none \
     --mincorr 0.5

With GMM clustering to find condition-specific relationships:

.. code-block:: bash

   kinc similarity \
     --input sample.emx.h5 \
     --cmx sample.cmx.h5 \
     --ccm sample.ccm.h5 \
     --corrmethod spearman \
     --clusmethod gmm \
     --minclus 1 \
     --maxclus 3 \
     --crit ICL \
     --mincorr 0.7

Step 4: Extract Network
~~~~~~~~~~~~~~~~~~~~~~~~

Extract edges above a correlation threshold:

.. code-block:: bash

   kinc extract \
     --input sample.cmx.h5 \
     --output sample_network.txt \
     --mincorr 0.8

The output file contains edges in tab-delimited format:

.. code-block:: text

   Source	Target	Cluster	Correlation
   Gene2	Gene1	0	0.952341
   Gene5	Gene2	0	0.887654
   Gene4	Gene3	0	0.823456

Using the Python API
--------------------

You can also use KINC-Python programmatically:

.. code-block:: python

   from kinc.analytics import import_emx, Similarity, extract_network

   # Import data
   import_emx("sample_data.txt", "sample.emx.h5")

   # Compute similarity
   sim = Similarity(
       input_file="sample.emx.h5",
       cmx_file="sample.cmx.h5",
       ccm_file="sample.ccm.h5",
       correlation_method='spearman',
       clustering_method='gmm',
       min_clusters=1,
       max_clusters=5,
       min_correlation=0.7
   )
   sim.run()

   # Extract network
   extract_network(
       cmx_file="sample.cmx.h5",
       output_file="network.txt",
       min_correlation=0.8
   )

Next Steps
----------

Now that you've run a basic workflow:

- Read the :doc:`usage` guide for detailed parameter explanations
- Check out :doc:`examples` for real-world workflows
- Learn about :doc:`data_formats` to understand KINC file structures
- See :doc:`api` for Python API documentation

Common Use Cases
----------------

Traditional Correlation Network
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a simple correlation network without clustering:

.. code-block:: bash

   kinc similarity \
     --input data.emx.h5 \
     --cmx data.cmx.h5 \
     --ccm data.ccm.h5 \
     --corrmethod pearson \
     --clusmethod none \
     --mincorr 0.85

Condition-Specific Network (GMM)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To find condition-specific correlations:

.. code-block:: bash

   kinc similarity \
     --input data.emx.h5 \
     --cmx data.cmx.h5 \
     --ccm data.ccm.h5 \
     --corrmethod spearman \
     --clusmethod gmm \
     --minclus 1 \
     --maxclus 5 \
     --crit ICL \
     --preout true \
     --postout true \
     --mincorr 0.85

Large Dataset on HPC
~~~~~~~~~~~~~~~~~~~~

For large datasets on HPC clusters with MPI:

.. code-block:: bash

   # Run with MPI (8 processes)
   mpirun -np 8 kinc similarity \
     --input large_data.emx.h5 \
     --cmx large_data.cmx.h5 \
     --ccm large_data.ccm.h5 \
     --corrmethod spearman \
     --clusmethod gmm \
     --mincorr 0.85

Tips and Best Practices
------------------------

**Choosing Correlation Method:**

- **Pearson**: For linear relationships, faster computation
- **Spearman**: For monotonic relationships, more robust to outliers

**Choosing Cluster Range:**

- Start with ``--minclus 1 --maxclus 5``
- Increase max if you expect many conditions
- Higher max increases computation time

**Setting Correlation Threshold:**

- Use ``--mincorr 0.85`` for stringent networks
- Use ``--mincorr 0.70`` for exploratory analysis
- Lower thresholds produce larger networks

**Memory and Performance:**

- For large datasets (>20K genes), use MPI or Dask
- Consider using ``--minsamp`` to require minimum samples
- Use ``--preout`` and ``--postout`` to remove outliers

Troubleshooting
---------------

**"Insufficient samples" warnings:**

Increase minimum samples: ``--minsamp 20``

**Very long computation time:**

- Use parallel processing (MPI or Dask)
- Reduce cluster range: ``--maxclus 3``
- Increase correlation threshold: ``--mincorr 0.9``

**Memory errors:**

- Process in chunks (MPI mode)
- Reduce number of genes in input
- Use a machine with more RAM

See :doc:`troubleshooting` for more help.
