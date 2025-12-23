Usage
=====

This page provides detailed information about using KINC-Python command-line tools and Python API.

Command-Line Interface
----------------------

KINC-Python provides a command-line tool ``kinc`` with several subcommands for different analytics.

General Syntax
~~~~~~~~~~~~~~

.. code-block:: bash

   kinc <command> [options]

Available Commands
~~~~~~~~~~~~~~~~~~

- ``import-emx``: Import gene expression matrix from text
- ``export-emx``: Export expression matrix to text
- ``similarity``: Compute pairwise gene correlations
- ``extract``: Extract network edges to file

Get help for any command:

.. code-block:: bash

   kinc <command> --help

import-emx
----------

Import a gene expression matrix from a text file to HDF5 format.

**Syntax:**

.. code-block:: bash

   kinc import-emx --input INPUT --output OUTPUT [options]

**Required Parameters:**

- ``--input``: Input text file (tab or space delimited)
- ``--output``: Output HDF5 file (.h5)

**Optional Parameters:**

- ``--nan``: Token for missing values (default: NA)
- ``--no-header``: Input file has no header row
- ``--no-rownames``: Input file has no row names (gene IDs)

**Examples:**

Import a standard GEM file:

.. code-block:: bash

   kinc import-emx --input data.txt --output data.emx.h5

Import with custom NA token:

.. code-block:: bash

   kinc import-emx --input data.txt --output data.emx.h5 --nan "N/A"

Import file without header:

.. code-block:: bash

   kinc import-emx --input data.txt --output data.emx.h5 --no-header

**Input File Format:**

The input file should be tab or space delimited with:

.. code-block:: text

   	Sample1	Sample2	Sample3	Sample4
   Gene1	1.23	2.45	1.78	2.01
   Gene2	3.45	4.12	3.89	4.23
   Gene3	0.56	0.78	0.67	0.82

export-emx
----------

Export an expression matrix from HDF5 to text format.

**Syntax:**

.. code-block:: bash

   kinc export-emx --input INPUT --output OUTPUT [options]

**Required Parameters:**

- ``--input``: Input HDF5 file (.h5)
- ``--output``: Output text file

**Optional Parameters:**

- ``--nan``: Token for missing values (default: NA)

**Example:**

.. code-block:: bash

   kinc export-emx --input data.emx.h5 --output data.txt

similarity
----------

Compute pairwise gene correlations with optional GMM clustering.

This is the core KINC analytic for network construction.

**Syntax:**

.. code-block:: bash

   kinc similarity --input INPUT --cmx CMX --ccm CCM [options]

**Required Parameters:**

- ``--input``: Input expression matrix (HDF5)
- ``--cmx``: Output correlation matrix (HDF5)
- ``--ccm``: Output cluster composition matrix (HDF5)

**Correlation Method:**

- ``--corrmethod``: Correlation method (default: pearson)

  - ``pearson``: Pearson correlation (linear relationships)
  - ``spearman``: Spearman rank correlation (monotonic relationships)

**Clustering Method:**

- ``--clusmethod``: Clustering method (default: none)

  - ``none``: No clustering, one correlation per pair
  - ``gmm``: Gaussian Mixture Model clustering

**GMM Parameters:**

- ``--minclus``: Minimum clusters to test (default: 1)
- ``--maxclus``: Maximum clusters to test (default: 5)
- ``--crit``: Model selection criterion (default: ICL)

  - ``AIC``: Akaike Information Criterion
  - ``BIC``: Bayesian Information Criterion
  - ``ICL``: Integrated Complete-data Likelihood

**Sample Filtering:**

- ``--minexpr``: Minimum expression threshold (default: -inf)
- ``--maxexpr``: Maximum expression threshold (default: +inf)
- ``--minsamp``: Minimum samples required (default: 30)

**Outlier Removal:**

- ``--preout``: Remove pre-clustering outliers (default: true)
- ``--postout``: Remove post-clustering outliers (default: true)

**Correlation Filtering:**

- ``--mincorr``: Minimum |correlation| to save (default: 0.5)
- ``--maxcorr``: Maximum |correlation| to save (default: 1.0)

**Performance:**

- ``--jobs``: Number of parallel jobs (default: 1)

**Examples:**

Basic Pearson correlation without clustering:

.. code-block:: bash

   kinc similarity \
     --input data.emx.h5 \
     --cmx data.cmx.h5 \
     --ccm data.ccm.h5 \
     --corrmethod pearson \
     --clusmethod none

GMM-based network with Spearman correlation:

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
     --mincorr 0.85

With custom expression thresholds:

.. code-block:: bash

   kinc similarity \
     --input data.emx.h5 \
     --cmx data.cmx.h5 \
     --ccm data.ccm.h5 \
     --corrmethod pearson \
     --minexpr 1.0 \
     --maxexpr 15.0 \
     --minsamp 25

extract
-------

Extract network edges from a correlation matrix to text file.

**Syntax:**

.. code-block:: bash

   kinc extract --input INPUT --output OUTPUT [options]

**Required Parameters:**

- ``--input``: Input correlation matrix (HDF5)
- ``--output``: Output network file

**Optional Parameters:**

- ``--ccm``: Cluster composition matrix (HDF5, optional)
- ``--mincorr``: Minimum |correlation| threshold (default: 0.0)
- ``--maxcorr``: Maximum |correlation| threshold (default: 1.0)
- ``--format``: Output format (default: text)

  - ``text``: Full format with sample masks
  - ``tidy``: Tidy format without sample masks
  - ``minimal``: Minimal format (only first cluster)

**Examples:**

Extract all edges above 0.85 correlation:

.. code-block:: bash

   kinc extract \
     --input data.cmx.h5 \
     --output network.txt \
     --mincorr 0.85

Extract with sample masks:

.. code-block:: bash

   kinc extract \
     --input data.cmx.h5 \
     --ccm data.ccm.h5 \
     --output network.txt \
     --mincorr 0.85 \
     --format text

Tidy format for easy import to R:

.. code-block:: bash

   kinc extract \
     --input data.cmx.h5 \
     --output network.txt \
     --format tidy

**Output Formats:**

*Text format (with CCM):*

.. code-block:: text

   Source	Target	Cluster	Correlation	SampleMask
   Gene1	Gene2	0	0.923	11111000000
   Gene1	Gene3	0	0.876	11110000000
   Gene2	Gene3	0	0.901	01111000000

*Tidy format:*

.. code-block:: text

   Gene1	Gene2	Cluster	Correlation
   Gene1	Gene2	0	0.923
   Gene1	Gene3	0	0.876
   Gene2	Gene3	0	0.901

*Minimal format (first cluster only):*

.. code-block:: text

   Gene1	Gene2	Correlation
   Gene1	Gene2	0.923
   Gene1	Gene3	0.876
   Gene2	Gene3	0.901

Python API
----------

KINC-Python can be used programmatically through its Python API.

Importing Modules
~~~~~~~~~~~~~~~~~

.. code-block:: python

   from kinc.analytics import import_emx, export_emx, Similarity, extract_network
   from kinc.data import ExpressionMatrix, CorrelationMatrix, CCMatrix

Importing Expression Matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import_emx(
       input_file="data.txt",
       output_file="data.emx.h5",
       nan_token="NA",
       has_header=True,
       has_row_names=True
   )

Computing Similarity
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   sim = Similarity(
       input_file="data.emx.h5",
       cmx_file="data.cmx.h5",
       ccm_file="data.ccm.h5",
       correlation_method='spearman',
       clustering_method='gmm',
       min_clusters=1,
       max_clusters=5,
       criterion='ICL',
       min_samples=30,
       preout=True,
       postout=True,
       min_correlation=0.85,
       max_correlation=1.0
   )

   # Run computation
   sim.run(n_jobs=1)

Extracting Network
~~~~~~~~~~~~~~~~~~

.. code-block:: python

   extract_network(
       cmx_file="data.cmx.h5",
       output_file="network.txt",
       ccm_file="data.ccm.h5",
       min_correlation=0.85,
       output_format='text'
   )

Accessing Data Directly
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from kinc.data import ExpressionMatrix, CorrelationMatrix

   # Read expression matrix
   with ExpressionMatrix("data.emx.h5", mode='r') as emx:
       print(f"Genes: {emx.n_genes}")
       print(f"Samples: {emx.n_samples}")

       # Read first gene
       expr = emx.read_expression(0)
       print(f"Gene 1 expression: {expr}")

   # Read correlation matrix
   with CorrelationMatrix("data.cmx.h5", mode='r') as cmx:
       # Iterate over all pairs
       for gene_x, gene_y, correlations in cmx.iter_pairs():
           print(f"{gene_x}, {gene_y}: {correlations}")

See :doc:`api` for complete API documentation.

Best Practices
--------------

**For Small Datasets (<1K genes):**

- Use serial mode (``--jobs 1``)
- Can run on laptop/workstation
- GMM clustering is fast enough

**For Medium Datasets (1K-10K genes):**

- Use multi-core parallelization (``--jobs 4`` or more)
- Consider reducing cluster range (``--maxclus 3``)
- May take several hours

**For Large Datasets (>10K genes):**

- Use MPI on HPC cluster
- Consider higher correlation threshold (``--mincorr 0.90``)
- May take days for very large datasets

**Choosing Model Selection Criterion:**

- **ICL** (default): Best for finding distinct clusters
- **BIC**: Good balance between fit and complexity
- **AIC**: May select more clusters

**Outlier Removal:**

- Enable both ``--preout`` and ``--postout`` for cleaner networks
- Disable if you want to preserve all data points
- Outliers are marked but not deleted from CCMatrix

See :doc:`performance` for optimization tips.
