
.. figure:: images/kinc.png
   :alt: KINC Logo

.. image:: https://zenodo.org/badge/71836133.svg
   :target: https://zenodo.org/badge/latestdoi/71836133

Welcome to the KINC documentation!
==================================

KINC is a C++ application for constructing a gene co-expression network (GCN) from a gene expression matrix (GEM). KINC implements several algorithms that are required for this process:

Pairwise Clustering

- Gaussian mixture models

Pairwise Correlation

- Pearson
- Spearman

Thresholding

- Power-law
- Random matrix theory

KINC is built with `ACE <https://github.com/SystemsGenetics/ACE>`__, a framework which provides mechanisms for large-scale heterogeneous computing and data management. As such, KINC can be run in a variety of compute configurations, including single-CPU / single-GPU and multi-CPU / multi-GPU, and KINC uses its own binary file formats to represent the data objects that it consumes and produces. Each of these binary formats can be exported to a plain-text format for use in other applications.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
