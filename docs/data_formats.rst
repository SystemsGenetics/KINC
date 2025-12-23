Data Formats
============

KINC-Python uses HDF5 format for all data storage, providing efficient, portable, and self-documenting files.

Why HDF5?
---------

HDF5 (Hierarchical Data Format version 5) offers several advantages:

- **Efficient**: Binary format with compression
- **Portable**: Works across platforms
- **Self-documenting**: Metadata embedded in file
- **Standard**: Widely used in scientific computing
- **Easy to inspect**: Tools like h5py, HDFView, h5dump
- **Interoperable**: Works with R, MATLAB, Julia, etc.

Expression Matrix (.emx.h5)
---------------------------

Stores gene expression data (genes × samples).

File Structure
~~~~~~~~~~~~~~

.. code-block:: text

   /
   ├── expressions: (n_genes, n_samples) float32 dataset
   ├── genes: (n_genes,) string dataset
   ├── samples: (n_samples,) string dataset
   └── metadata: group
       ├── version: "3.4.2" (attribute)
       ├── n_genes: integer (attribute)
       └── n_samples: integer (attribute)

Accessing in Python
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import h5py

   with h5py.File('data.emx.h5', 'r') as f:
       # Read metadata
       n_genes = f['metadata'].attrs['n_genes']
       n_samples = f['metadata'].attrs['n_samples']

       # Read gene names
       genes = f['genes'][:]

       # Read expression data
       expressions = f['expressions'][:]  # Full matrix
       gene_expr = f['expressions'][0, :]  # First gene

Using KINC API:

.. code-block:: python

   from kinc.data import ExpressionMatrix

   with ExpressionMatrix('data.emx.h5', mode='r') as emx:
       print(f"Genes: {emx.n_genes}")
       print(f"Samples: {emx.n_samples}")
       expr = emx.read_expression(0)  # First gene

Correlation Matrix (.cmx.h5)
----------------------------

Stores pairwise gene correlations with multiple clusters per pair.

File Structure
~~~~~~~~~~~~~~

.. code-block:: text

   /
   ├── pairs: group
   │   ├── pair_00000000: group
   │   │   ├── indices: (2,) int32 [gene_x, gene_y]
   │   │   └── correlations: (n_clusters,) float32
   │   ├── pair_00000001: group
   │   │   ├── indices: (2,) int32
   │   │   └── correlations: (n_clusters,) float32
   │   └── ...
   ├── genes: (n_genes,) string dataset
   └── metadata: group
       ├── version: "3.4.2"
       ├── n_genes: integer
       ├── correlation_method: "pearson" or "spearman"
       └── max_clusters: integer

Sparse Storage
~~~~~~~~~~~~~~

Only gene pairs with valid correlations are stored (not all N×N pairs).

Pair IDs are calculated as:

.. code-block:: python

   pair_id = gene_x * (gene_x - 1) // 2 + gene_y

where ``gene_x > gene_y`` (lower triangle only).

Accessing in Python
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from kinc.data import CorrelationMatrix

   with CorrelationMatrix('data.cmx.h5', mode='r') as cmx:
       # Read specific pair
       corrs = cmx.read_pair(gene_x=5, gene_y=2)

       # Iterate all pairs
       for gene_x, gene_y, correlations in cmx.iter_pairs():
           print(f"Genes {gene_x}, {gene_y}: {correlations}")

CCMatrix (.ccm.h5)
------------------

Stores cluster sample masks (which samples belong to each cluster).

File Structure
~~~~~~~~~~~~~~

.. code-block:: text

   /
   ├── pairs: group
   │   ├── pair_00000000: group
   │   │   ├── indices: (2,) int32
   │   │   └── sample_masks: (n_clusters, n_samples) int8
   │   └── ...
   ├── genes: (n_genes,) string dataset
   ├── samples: (n_samples,) string dataset
   └── metadata: group
       ├── version: "3.4.2"
       ├── n_genes: integer
       └── n_samples: integer

Sample Mask Values
~~~~~~~~~~~~~~~~~~

Each sample has an integer label:

- ``0-15``: Cluster assignment (0 = cluster 0, 1 = cluster 1, etc.)
- ``-6``: Below/above expression threshold
- ``-7``: Pre-clustering outlier
- ``-8``: Post-clustering outlier
- ``-9``: Missing value (NaN)

Accessing in Python
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from kinc.data import CCMatrix

   with CCMatrix('data.ccm.h5', mode='r') as ccm:
       # Read sample masks for a pair
       masks = ccm.read_pair(gene_x=5, gene_y=2)

       # masks shape: (n_clusters, n_samples)
       # masks[0, :] = sample mask for cluster 0
       # masks[1, :] = sample mask for cluster 1

Inspecting HDF5 Files
----------------------

Command-Line Tools
~~~~~~~~~~~~~~~~~~

List file structure:

.. code-block:: bash

   h5ls -r data.emx.h5

Dump file contents:

.. code-block:: bash

   h5dump data.emx.h5

View with HDFView
~~~~~~~~~~~~~~~~~

Download `HDFView <https://www.hdfgroup.org/downloads/hdfview/>`_ for a graphical interface.

Using Python h5py
~~~~~~~~~~~~~~~~~

.. code-block:: python

   import h5py

   with h5py.File('data.emx.h5', 'r') as f:
       # List all datasets
       def print_structure(name, obj):
           print(name)
       f.visititems(print_structure)

       # Read attributes
       for key, val in f['metadata'].attrs.items():
           print(f"{key}: {val}")

Converting to Other Formats
----------------------------

To CSV/TSV
~~~~~~~~~~

.. code-block:: python

   from kinc.analytics import export_emx

   export_emx('data.emx.h5', 'data.txt')

To NumPy Arrays
~~~~~~~~~~~~~~~

.. code-block:: python

   from kinc.data import ExpressionMatrix
   import numpy as np

   emx = ExpressionMatrix('data.emx.h5', mode='r')
   data = emx.read_expressions()  # NumPy array

   # Save as .npy
   np.save('data.npy', data)

To Pandas DataFrame
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from kinc.data import ExpressionMatrix
   import pandas as pd

   emx = ExpressionMatrix('data.emx.h5', mode='r')
   data = emx.read_expressions()

   df = pd.DataFrame(
       data,
       index=emx.gene_names,
       columns=emx.sample_names
   )

   # Save as CSV
   df.to_csv('data.csv')

   # Save as Parquet
   df.to_parquet('data.parquet')

File Size Comparison
--------------------

Typical file sizes for a 20,000 gene × 100 sample dataset:

.. list-table::
   :header-rows: 1
   :widths: 40 20 40

   * - Format
     - Size
     - Notes
   * - Text (uncompressed)
     - ~400 MB
     - Large, slow to read
   * - HDF5 (compressed)
     - ~80 MB
     - 5x smaller, fast random access
   * - NumPy .npy
     - ~160 MB
     - Fast but no metadata
   * - Parquet
     - ~100 MB
     - Good compression, pandas-friendly

For correlation matrices (sparse storage):

.. list-table::
   :header-rows: 1
   :widths: 40 20 40

   * - Pairs Stored
     - Size
     - Notes
   * - 1 million pairs
     - ~40 MB
     - Average 5 clusters/pair
   * - 10 million pairs
     - ~400 MB
     - Typical for 5K genes
   * - 100 million pairs
     - ~4 GB
     - Typical for 15K genes

Best Practices
--------------

**File Naming:**

- Use descriptive names: ``rice_drought_stress.emx.h5``
- Include dataset identifier
- Use consistent suffixes: ``.emx.h5``, ``.cmx.h5``, ``.ccm.h5``

**Storage:**

- Keep raw text files archived
- HDF5 files are the working format
- Back up correlation matrices (expensive to recompute)

**Compression:**

- KINC-Python uses gzip compression by default
- Balance: compression level 4 is good trade-off
- Higher levels save space but slow I/O

**Memory:**

- HDF5 allows chunked reading (doesn't load all in memory)
- Process data in chunks for large datasets
- Use iterators when possible

Compatibility
-------------

**With C++ KINC:**

KINC-Python uses HDF5, while C++ KINC uses custom binary formats.
They are **not directly compatible**, but you can:

1. Export from C++ KINC to text
2. Import to Python KINC from text

**With R:**

HDF5 files work seamlessly with R:

.. code-block:: r

   library(rhdf5)

   # Read expression matrix
   emx <- h5read("data.emx.h5", "expressions")
   genes <- h5read("data.emx.h5", "genes")
   samples <- h5read("data.emx.h5", "samples")

**With Other Tools:**

Most scientific software supports HDF5:

- MATLAB: ``h5read()``, ``h5write()``
- Julia: HDF5.jl package
- Java: HDF5 Java library
- Command-line: h5dump, h5ls, h5copy
