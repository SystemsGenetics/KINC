Troubleshooting
===============

There are many issues that can occur when using KINC, ranging from data management to software and hardware configuration. The following sections will attempt to address the most common issues.

Thresholding Large GEMs
~~~~~~~~~~~~~~~~~~~~~~~

There are several issues that emerge when applying KINC to large input GEMs. Aside from the large computational cost in the similarity step, the main issue with large GEMs is the thresholding step. The ``rmt`` analytic finds a suitable threshold by performing a chi-squared test for each threshold to determine whether the thresholded network is random or non-random. In particular, the chi-squared value should transition from a low value (< 100) to a high value (> 200), and the threshold at which this transition occurs is selected as the final threshold. However, if the similarity matrix contains many edges even at high correlations (> 0.95), it will likely cause the chi-squared value to be too large for this transition to ever occur, and a suitable threshold will not be found.

There are many approaches to dealing with this issue. In general, the best thing to do first is to visualize some properties of the data: the sample distribution of the GEM (consult the GEMprep tool), the correlation distribution of the similarity matrix (using ``visualize.py``), and pairwise scatter plots of the similarity matrix (using ``visualize.py``). These visualizations can help you to identify irregularities or noise in the data, from which you can decide how to remove them. For example, if the GEM contains many negative-valued expressions that are noisy, you can use the ``--minexpr`` option in the ``similarity`` analytic to exclude these expressions when computing the pairwise correlations, which may remove many high but meaningless correlations.
