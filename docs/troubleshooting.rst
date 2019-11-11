Troubleshooting
===============

There are many issues that can occur when using KINC, ranging from data management to software and hardware configuration. The following sections will attempt to address the most common issues.

Thresholding Large GEMs
~~~~~~~~~~~~~~~~~~~~~~~

There are several issues that emerge when applying KINC to large input GEMs. Aside from the large computational cost in the similarity step, the main issue with large GEMs is the thresholding step. The ``rmt`` analytic finds a suitable threshold by performing a chi-squared test for each threshold to determine whether the thresholded network is random or non-random. In particular, the chi-squared value should transition from a low value (< 100) to a high value (> 200), and the threshold at which this transition occurs is selected as the final threshold. However, if the similarity matrix contains many edges even at high correlations (> 0.95), it will likely cause the chi-squared value to be too large for this transition to ever occur, and a suitable threshold will not be found.

There are many approaches to dealing with this issue. In general, the best thing to do first is to visualize some properties of the data: the sample distribution of the GEM (consult the GEMprep tool), the correlation distribution of the similarity matrix (using ``visualize.py``), and pairwise scatter plots of the similarity matrix (using ``visualize.py``). These visualizations can help you to identify irregularities or noise in the data, from which you can decide how to remove them. For example, if the GEM contains many negative-valued expressions that are noisy, you can use the ``--minexpr`` option in the ``similarity`` analytic to exclude these expressions when computing the pairwise correlations, which may remove many high but meaningless correlations.

Performance Considerations
~~~~~~~~~~~~~~~~~~~~~~~~~~

Since KINC can be run with a variety of hardware configurations, including single-CPU, multi-CPU, single-GPU, and multi-GPU, there are several settings that control how KINC uses this hardware. In particular, the multi-GPU configuration for ``similarity`` is the most complex and uses all of the execution parameters. Here we describe each execution parameter and provide recommendations based on performance benchmarking and experience.

- **CUDA/OpenCL Thread Size**: Determines the number of worker threads per GPU. Increasing this value can increase performance by utilizing the GPU more fully, but setting this value too high can also decrease performance due to the overhead of switching between many threads. A safe value for this parameter is 2 threads, however on newer hardware it may be possible to use more threads and achieve better performance. This parameter is set using the ``threads`` option in the KINC settings.

- **MPI Work Block Size**: Determines the number of work items per MPI work block. It is effectively the maximum number of work items that a worker thread can process in parallel. In practice, the work block size does not affect performance so long as it is greater than or equal to the global work size, so the default value of 32,768 should work well. This parameter is set using the ``--bsize`` option in the ``similarity`` analytic.

- **Global Work Size**: Determines the number of work items that a worker thread processes in parallel on the GPU. It should be large enough to fully utilize the GPU, but setting it too large can also decrease performance due to global memory congestion and work imbalance on the GPU. In practice, the default value of 4096 seems to work the best. This parameter is set using the ``--gisze`` option in the ``similarity`` analytic.

- **Local Work Size**: Determines the OpenCL local work size (CUDA block size) of each GPU kernel. In general, the optimal value for this parameter depends heavily on the particular GPU kernel, but since all of the GPU kernels in KINC are memory-intensive, the local work size should be small to prevent global memory congestion. In practice, a value of 16 or 32 (the default) works the best. This parameter is set using the ``--lsize`` option in the ``similarity`` analytic.
