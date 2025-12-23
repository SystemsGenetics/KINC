Performance
===========

This page provides guidance on optimizing KINC-Python performance for different scenarios.

Performance Overview
--------------------

Typical performance for a 20,000 gene × 100 sample dataset:

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 30

   * - Mode
     - Time
     - Memory
     - Notes
   * - Serial (no GMM)
     - ~4 hours
     - 8 GB
     - Single core
   * - Serial (GMM)
     - ~12 hours
     - 10 GB
     - Single core
   * - Parallel 8 cores
     - ~2 hours
     - 12 GB
     - When implemented
   * - MPI 20 nodes
     - ~30 min
     - 40 GB total
     - When implemented
   * - GPU (CuPy)
     - ~1 hour
     - 16 GB GPU
     - When implemented

Factors Affecting Performance
------------------------------

Dataset Size
~~~~~~~~~~~~

Performance scales with:

- **Number of genes**: O(N²) pairs to compute
- **Number of samples**: Linear impact on correlation computation
- **Cluster complexity**: More clusters = slower GMM fitting

.. list-table::
   :header-rows: 1

   * - Genes
     - Pairs
     - Serial Time (est.)
   * - 1,000
     - 500K
     - ~10 minutes
   * - 5,000
     - 12.5M
     - ~2 hours
   * - 10,000
     - 50M
     - ~10 hours
   * - 20,000
     - 200M
     - ~48 hours

Algorithm Choices
~~~~~~~~~~~~~~~~~

**Correlation Method:**

- **Pearson**: Faster (~2x) - O(N) complexity
- **Spearman**: Slower - requires ranking O(N log N)

**Clustering Method:**

- **None**: Fast - single correlation per pair
- **GMM**: Slow - tests K=1 to maxclus, runs EM algorithm

**Model Selection:**

- **AIC**: Fastest - simple calculation
- **BIC**: Fast - similar to AIC
- **ICL**: Slower - requires entropy calculation

Optimization Strategies
-----------------------

For Small Datasets (<1K genes)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use serial mode with default settings:

.. code-block:: bash

   kinc similarity \
     --input data.emx.h5 \
     --cmx data.cmx.h5 \
     --ccm data.ccm.h5 \
     --corrmethod pearson \
     --clusmethod gmm \
     --mincorr 0.85

**Expected time**: Minutes to 1 hour

For Medium Datasets (1K-10K genes)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use multi-core parallelization (when implemented):

.. code-block:: bash

   kinc similarity \
     --input data.emx.h5 \
     --cmx data.cmx.h5 \
     --ccm data.ccm.h5 \
     --corrmethod pearson \
     --clusmethod gmm \
     --maxclus 3 \
     --mincorr 0.90 \
     --jobs 8

**Optimizations**:

- Reduce max clusters: ``--maxclus 3``
- Increase correlation threshold: ``--mincorr 0.90``
- Use multiple cores: ``--jobs 8``

**Expected time**: 1-4 hours

For Large Datasets (>10K genes)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use MPI on HPC cluster:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --nodes=5
   #SBATCH --ntasks=20
   #SBATCH --time=8:00:00

   module load kinc-python
   module load openmpi

   mpirun -np 20 kinc similarity \
     --input data.emx.h5 \
     --cmx data.cmx.h5 \
     --ccm data.ccm.h5 \
     --corrmethod spearman \
     --clusmethod gmm \
     --maxclus 3 \
     --mincorr 0.85

**Optimizations**:

- Use MPI across nodes
- Conservative cluster range
- High correlation threshold
- Consider preprocessing to filter genes

**Expected time**: 2-8 hours

Parameter Tuning
----------------

Correlation Threshold
~~~~~~~~~~~~~~~~~~~~~

Higher thresholds reduce computation and storage:

.. list-table::
   :header-rows: 1

   * - Threshold
     - Typical Edges
     - Savings
   * - 0.50
     - ~10M
     - Baseline
   * - 0.70
     - ~5M
     - 50%
   * - 0.85
     - ~1M
     - 90%
   * - 0.90
     - ~500K
     - 95%

**Recommendation**: Start with 0.85 for most analyses.

Cluster Range
~~~~~~~~~~~~~

Fewer clusters = faster computation:

.. list-table::
   :header-rows: 1

   * - Max Clusters
     - Models Tested
     - Time Impact
   * - 2
     - 2
     - ~2x faster
   * - 3
     - 3
     - ~1.5x faster
   * - 5 (default)
     - 5
     - Baseline
   * - 7
     - 7
     - ~1.4x slower

**Recommendation**: Use ``--maxclus 3`` unless you expect many conditions.

Minimum Samples
~~~~~~~~~~~~~~~

Higher minimum = fewer pairs computed:

.. code-block:: bash

   # Require at least 40 samples (stricter)
   kinc similarity ... --minsamp 40

Trade-off: Fewer false positives vs. more pairs excluded.

Memory Optimization
-------------------

Reduce Memory Usage
~~~~~~~~~~~~~~~~~~~

1. **Higher correlation threshold**:

.. code-block:: bash

   --mincorr 0.90  # Saves fewer pairs

2. **Fewer clusters**:

.. code-block:: bash

   --maxclus 3

3. **Close other applications**

4. **Use swap space** (Linux):

.. code-block:: bash

   sudo swapon /swapfile

Monitor Memory
~~~~~~~~~~~~~~

.. code-block:: bash

   # During execution
   watch -n 1 free -h

   # Or with Python
   import psutil
   print(f"Memory: {psutil.virtual_memory().percent}%")

Disk I/O Optimization
---------------------

Use Fast Storage
~~~~~~~~~~~~~~~~

- **Best**: NVMe SSD
- **Good**: SATA SSD
- **Slow**: HDD, Network file systems

For HPC, use local scratch space:

.. code-block:: bash

   # Copy to local scratch
   cp data.emx.h5 $TMPDIR/
   cd $TMPDIR

   # Run KINC
   kinc similarity ...

   # Copy results back
   cp *.cmx.h5 *.ccm.h5 $HOME/results/

HDF5 Compression
~~~~~~~~~~~~~~~~

KINC-Python uses gzip compression (level 4) by default. This balances:

- **Compression**: ~5x reduction in file size
- **Speed**: Minimal overhead for reading/writing

Parallel Processing
-------------------

Multi-Core (Future)
~~~~~~~~~~~~~~~~~~~

When implemented, use:

.. code-block:: bash

   kinc similarity ... --jobs $(nproc)

Or specify cores:

.. code-block:: bash

   kinc similarity ... --jobs 8

MPI (Future)
~~~~~~~~~~~~

For HPC clusters:

.. code-block:: bash

   #!/bin/bash
   #SBATCH --nodes=5
   #SBATCH --ntasks-per-node=4

   mpirun -np 20 kinc similarity ...

**Best practices**:

- Use 4-8 tasks per node
- Match tasks to available cores
- Ensure sufficient memory per task

GPU Acceleration (Future)
~~~~~~~~~~~~~~~~~~~~~~~~~~

When CuPy support is added:

.. code-block:: bash

   kinc similarity ... --gpu

**Requirements**:

- NVIDIA GPU with CUDA
- CuPy installed
- Sufficient GPU memory

Benchmarking
------------

Measure Performance
~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Time execution
   time kinc similarity ...

   # Or with more detail
   /usr/bin/time -v kinc similarity ...

Create Test Dataset
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import numpy as np
   import pandas as pd

   # Small test (100 genes)
   n_genes = 100
   n_samples = 50

   data = np.random.randn(n_genes, n_samples)
   df = pd.DataFrame(
       data,
       index=[f"Gene_{i}" for i in range(n_genes)],
       columns=[f"Sample_{i}" for i in range(n_samples)]
   )
   df.to_csv("test_small.txt", sep='\t')

Run benchmark:

.. code-block:: bash

   kinc import-emx --input test_small.txt --output test.emx.h5
   time kinc similarity --input test.emx.h5 --cmx test.cmx.h5 --ccm test.ccm.h5

Profiling
---------

Profile Python Code
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import cProfile
   import pstats
   from kinc.analytics import Similarity

   # Profile similarity computation
   profiler = cProfile.Profile()
   profiler.enable()

   sim = Similarity("data.emx.h5", "data.cmx.h5", "data.ccm.h5")
   sim.run()

   profiler.disable()

   # Print stats
   stats = pstats.Stats(profiler)
   stats.sort_stats('cumulative')
   stats.print_stats(20)

Identify Bottlenecks
~~~~~~~~~~~~~~~~~~~~

Common bottlenecks:

1. **GMM fitting**: Most time-consuming
2. **Correlation computation**: Moderate
3. **File I/O**: Usually fast with HDF5

Comparison with C++ KINC
-------------------------

.. list-table::
   :header-rows: 1

   * - Feature
     - C++ KINC
     - Python KINC
     - Notes
   * - Serial Performance
     - Baseline
     - 2-5x slower
     - NumPy overhead
   * - Memory Usage
     - Lower
     - Similar
     - HDF5 is efficient
   * - GPU (CUDA)
     - Very fast
     - Pending
     - Custom vs CuPy
   * - MPI
     - Very fast
     - Pending
     - Native vs mpi4py
   * - Ease of Use
     - Moderate
     - Easy
     - Python simpler

**Conclusion**: Python version is slower but more accessible. For very large datasets, consider C++ KINC or wait for parallel Python support.

Best Practices Summary
----------------------

1. **Start small**: Test workflow on small dataset first
2. **Choose appropriate method**: Pearson for speed, Spearman for robustness
3. **Use conservative settings**: ``--maxclus 3 --mincorr 0.85``
4. **Monitor resources**: Check memory and CPU usage
5. **Use appropriate hardware**: SSDs for I/O, sufficient RAM
6. **Plan for scale**: Use MPI for large datasets (when available)
7. **Filter early**: Preprocess to remove low-variance genes
8. **Save intermediate results**: Keep expression matrix, can rerun similarity

Performance Checklist
----------------------

Before running on large dataset:

- [ ] Test on small subset (100-1000 genes)
- [ ] Verify settings work correctly
- [ ] Estimate time based on test run
- [ ] Ensure sufficient disk space (2-3x dataset size)
- [ ] Check available memory (8+ GB recommended)
- [ ] Plan for interruptions (use screen/tmux)
- [ ] Back up input data
- [ ] Monitor first hour of execution

See Also
--------

- :doc:`troubleshooting` for common performance issues
- :doc:`usage` for parameter details
- :doc:`examples` for real-world workflows
