Troubleshooting
===============

This page provides solutions to common problems when using KINC-Python.

Installation Issues
-------------------

Package Not Found
~~~~~~~~~~~~~~~~~

**Error**: ``Could not find a version that satisfies the requirement kinc-python``

**Solution**:

.. code-block:: bash

   # Update pip
   pip install --upgrade pip

   # Install from GitHub
   pip install git+https://github.com/SystemsGenetics/KINC-Python.git

HDF5 Installation Fails
~~~~~~~~~~~~~~~~~~~~~~~~

**Error**: ``h5py installation failed`` or ``Could not find HDF5``

**Solution**:

Install system HDF5 library first:

.. code-block:: bash

   # Ubuntu/Debian
   sudo apt install libhdf5-dev

   # macOS
   brew install hdf5

   # Then reinstall h5py
   pip install --no-cache-dir h5py

CuPy Installation Issues
~~~~~~~~~~~~~~~~~~~~~~~~~

**Error**: ``cupy installation failed``

**Solution**:

1. Verify CUDA is installed:

.. code-block:: bash

   nvcc --version
   nvidia-smi

2. Install matching CuPy version:

.. code-block:: bash

   # For CUDA 12.x
   pip install cupy-cuda12x

   # For CUDA 11.x
   pip install cupy-cuda11x

3. See `CuPy installation guide <https://docs.cupy.dev/en/stable/install.html>`_

Runtime Issues
--------------

Import Errors
~~~~~~~~~~~~~

**Error**: ``No module named 'kinc'``

**Solutions**:

1. Check installation:

.. code-block:: bash

   pip list | grep kinc

2. Verify Python environment:

.. code-block:: bash

   which python
   which pip

3. Reinstall in correct environment:

.. code-block:: bash

   pip install -e .

Memory Errors
~~~~~~~~~~~~~

**Error**: ``MemoryError`` or ``Out of memory``

**Solutions**:

1. **Reduce batch size** (when implemented):

.. code-block:: bash

   kinc similarity --input data.emx.h5 ... --bsize 10000

2. **Use MPI** to distribute across nodes

3. **Increase correlation threshold**:

.. code-block:: bash

   kinc similarity ... --mincorr 0.90  # Saves fewer pairs

4. **Close other applications**

5. **Use swap space** (Linux):

.. code-block:: bash

   sudo fallocate -l 16G /swapfile
   sudo chmod 600 /swapfile
   sudo mkswap /swapfile
   sudo swapon /swapfile

File Format Issues
------------------

HDF5 File Corrupted
~~~~~~~~~~~~~~~~~~~

**Error**: ``Unable to open file`` or ``File signature not found``

**Solutions**:

1. **Check file integrity**:

.. code-block:: bash

   h5dump data.emx.h5 > /dev/null

2. **Re-import from text** if possible

3. **Check disk space**:

.. code-block:: bash

   df -h

Invalid Expression Values
~~~~~~~~~~~~~~~~~~~~~~~~~

**Error**: ``Invalid values in expression matrix``

**Solutions**:

1. **Check for inf/nan**:

.. code-block:: python

   import numpy as np
   from kinc.data import ExpressionMatrix

   with ExpressionMatrix("data.emx.h5", mode='r') as emx:
       data = emx.read_expressions()
       print(f"NaN count: {np.isnan(data).sum()}")
       print(f"Inf count: {np.isinf(data).sum()}")

2. **Clean data before import**:

.. code-block:: python

   import pandas as pd

   df = pd.read_csv("data.txt", sep='\t', index_col=0)
   df = df.replace([np.inf, -np.inf], np.nan)
   df = df.fillna(0)  # Or drop: df = df.dropna()
   df.to_csv("data_clean.txt", sep='\t')

Computation Issues
------------------

Very Slow Performance
~~~~~~~~~~~~~~~~~~~~~

**Problem**: Similarity computation taking too long

**Solutions**:

1. **Use parallel processing**:

.. code-block:: bash

   # Multi-core (when implemented)
   kinc similarity ... --jobs 8

   # MPI on cluster
   mpirun -np 20 kinc similarity ...

2. **Reduce cluster range**:

.. code-block:: bash

   kinc similarity ... --maxclus 3  # Instead of 5

3. **Increase correlation threshold**:

.. code-block:: bash

   kinc similarity ... --mincorr 0.90

4. **Use Pearson instead of Spearman**:

.. code-block:: bash

   kinc similarity ... --corrmethod pearson  # Faster than spearman

5. **Reduce number of genes**:

Filter to most variable genes before importing.

Insufficient Samples Warning
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Warning**: ``Pair skipped: insufficient samples``

**Causes**:
- Too many NaN values
- Expression thresholds too stringent
- Not enough valid samples after outlier removal

**Solutions**:

1. **Reduce minimum samples**:

.. code-block:: bash

   kinc similarity ... --minsamp 20  # Instead of 30

2. **Relax expression thresholds**:

.. code-block:: bash

   kinc similarity ... --minexpr -10 --maxexpr 20

3. **Disable outlier removal**:

.. code-block:: bash

   kinc similarity ... --preout false --postout false

4. **Check data quality**:

.. code-block:: python

   from kinc.data import ExpressionMatrix
   import numpy as np

   emx = ExpressionMatrix("data.emx.h5", mode='r')
   data = emx.read_expressions()

   # Count valid samples per gene
   valid_per_gene = np.sum(~np.isnan(data), axis=1)
   print(f"Min valid samples: {valid_per_gene.min()}")
   print(f"Genes with <30 samples: {np.sum(valid_per_gene < 30)}")

GMM Clustering Failures
~~~~~~~~~~~~~~~~~~~~~~~~

**Problem**: GMM fails to converge or selects wrong K

**Solutions**:

1. **Try different criterion**:

.. code-block:: bash

   # BIC is more conservative than AIC
   kinc similarity ... --crit BIC

2. **Adjust cluster range**:

.. code-block:: bash

   # More conservative
   kinc similarity ... --minclus 1 --maxclus 3

3. **Check data distribution**:

GMM works best when samples form distinct groups.

Network Extraction Issues
--------------------------

No Edges Extracted
~~~~~~~~~~~~~~~~~~

**Problem**: Extract produces empty file

**Causes**:
- Correlation threshold too high
- No pairs meet criteria

**Solutions**:

1. **Lower threshold**:

.. code-block:: bash

   kinc extract ... --mincorr 0.70  # Instead of 0.90

2. **Check CMX file**:

.. code-block:: python

   from kinc.data import CorrelationMatrix

   with CorrelationMatrix("data.cmx.h5", mode='r') as cmx:
       print(f"Total pairs: {cmx.pair_count}")

3. **Inspect correlations**:

.. code-block:: python

   import numpy as np

   all_corrs = []
   for _, _, correlations in cmx.iter_pairs():
       all_corrs.extend(correlations[~np.isnan(correlations)])

   print(f"Min corr: {np.min(all_corrs)}")
   print(f"Max corr: {np.max(all_corrs)}")
   print(f"Median: {np.median(all_corrs)}")

Python API Issues
-----------------

Context Manager Errors
~~~~~~~~~~~~~~~~~~~~~~~

**Error**: ``File not opened``

**Solution**: Use ``with`` statement:

.. code-block:: python

   # Wrong
   emx = ExpressionMatrix("data.emx.h5", mode='r')
   data = emx.read_expressions()  # Error!

   # Correct
   with ExpressionMatrix("data.emx.h5", mode='r') as emx:
       data = emx.read_expressions()

Type Errors
~~~~~~~~~~~

**Error**: ``TypeError: expected str, got int``

**Solution**: Convert types explicitly:

.. code-block:: python

   # Wrong
   gene_idx = "0"

   # Correct
   gene_idx = 0

HPC/Cluster Issues
------------------

MPI Not Found
~~~~~~~~~~~~~

**Error**: ``mpirun: command not found``

**Solution**:

Load MPI module:

.. code-block:: bash

   module load openmpi
   mpirun -np 4 kinc similarity ...

Permission Denied
~~~~~~~~~~~~~~~~~

**Error**: ``Permission denied`` on cluster

**Solutions**:

1. **Check file permissions**:

.. code-block:: bash

   ls -la data.emx.h5

2. **Use scratch space**:

.. code-block:: bash

   # Copy to scratch
   cp data.emx.h5 $SCRATCH/
   cd $SCRATCH
   kinc similarity ...

3. **Check quota**:

.. code-block:: bash

   quota -s

Module Load Errors
~~~~~~~~~~~~~~~~~~

**Error**: ``Module not found``

**Solution**:

.. code-block:: bash

   # List available modules
   module avail

   # Load required modules
   module load python/3.10
   module load gcc/11.2

Performance Optimization
------------------------

Speed Up Computation
~~~~~~~~~~~~~~~~~~~~

1. **Use appropriate correlation method**:

   - Pearson: ~2x faster than Spearman
   - Spearman: More robust to outliers

2. **Optimize cluster settings**:

   - Start with ``--maxclus 3``
   - Increase only if needed

3. **Filter genes beforehand**:

   - Keep most variable genes only
   - Remove low-expression genes

4. **Use appropriate thresholds**:

   - Higher ``--mincorr`` = fewer pairs saved
   - Reduces memory and disk usage

Reduce Memory Usage
~~~~~~~~~~~~~~~~~~~

1. **Process in chunks** (when MPI implemented)

2. **Increase correlation threshold**

3. **Use fewer clusters**: ``--maxclus 3``

4. **Close unused applications**

Getting Help
------------

If you encounter issues not covered here:

1. **Check GitHub Issues**:

   https://github.com/SystemsGenetics/KINC-Python/issues

2. **Search for similar problems**

3. **Create new issue** with:

   - KINC-Python version
   - Python version
   - Operating system
   - Full error message
   - Minimal example to reproduce

4. **Contact developers**:

   stephen.ficklin@wsu.edu

Debugging Tips
--------------

Enable Verbose Output
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import logging
   logging.basicConfig(level=logging.DEBUG)

Check File Contents
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # List HDF5 structure
   h5ls -r data.emx.h5

   # Dump contents
   h5dump data.emx.h5 | head -100

Validate Data
~~~~~~~~~~~~~

.. code-block:: python

   from kinc.data import ExpressionMatrix
   import numpy as np

   with ExpressionMatrix("data.emx.h5", mode='r') as emx:
       data = emx.read_expressions()

       print(f"Shape: {data.shape}")
       print(f"Dtype: {data.dtype}")
       print(f"NaN count: {np.isnan(data).sum()}")
       print(f"Range: [{data.min()}, {data.max()}]")

Test with Small Dataset
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Create small test dataset
   import numpy as np
   import pandas as pd

   data = np.random.randn(10, 20)
   df = pd.DataFrame(
       data,
       index=[f"Gene_{i}" for i in range(10)],
       columns=[f"Sample_{i}" for i in range(20)]
   )
   df.to_csv("test.txt", sep='\t')

   # Test workflow
   from kinc.analytics import import_emx, Similarity

   import_emx("test.txt", "test.emx.h5")
   sim = Similarity("test.emx.h5", "test.cmx.h5", "test.ccm.h5")
   sim.run()
