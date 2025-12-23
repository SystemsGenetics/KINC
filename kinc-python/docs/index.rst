.. image:: https://zenodo.org/badge/71836133.svg
   :target: https://zenodo.org/badge/latestdoi/71836133

Welcome to KINC-Python Documentation!
======================================

The Knowledge Independent Network Construction (KINC) Python implementation is a complete rewrite of the original C++ KINC software for constructing `Gene Co-expression Networks (GCN) <https://en.wikipedia.org/wiki/Gene_co-expression_network>`_ from gene expression data.

KINC-Python leverages modern Python scientific libraries to provide:

- **GPU-accelerated** pairwise correlation computation (Pearson, Spearman)
- **Gaussian Mixture Model (GMM)** clustering for condition-specific networks
- **MPI support** for distributed computing on HPC clusters
- **HDF5 data formats** for efficient storage and interoperability
- **Command-line and Python API** for flexible workflows

.. toctree::
   :maxdepth: 3
   :caption: Contents:

   about
   installation
   quickstart
   usage
   api
   examples
   data_formats
   performance
   troubleshooting
   cite
   contributing

Key Features
------------

- **Modern Python Stack**: Built on NumPy, SciPy, scikit-learn, and HDF5
- **Statistically Equivalent**: Results match C++ KINC within numerical tolerance
- **Easy to Use**: Simple command-line interface and Python API
- **Extensible**: Easy to add new correlation methods or clustering algorithms
- **Well Documented**: Comprehensive documentation and examples

Quick Example
-------------

.. code-block:: bash

   # Import gene expression matrix
   kinc import-emx --input data.txt --output data.emx.h5

   # Compute similarity with GMM clustering
   kinc similarity \
     --input data.emx.h5 \
     --cmx data.cmx.h5 \
     --ccm data.ccm.h5 \
     --corrmethod spearman \
     --clusmethod gmm \
     --minclus 1 \
     --maxclus 5

   # Extract network edges
   kinc extract --input data.cmx.h5 --output network.txt --mincorr 0.85

Installation
------------

.. code-block:: bash

   pip install kinc-python

Or for development:

.. code-block:: bash

   git clone https://github.com/SystemsGenetics/KINC-Python.git
   cd KINC-Python
   pip install -e .

For GPU support:

.. code-block:: bash

   pip install kinc-python[gpu]

See :doc:`installation` for detailed instructions.

Citation
--------

If you use KINC-Python in your research, please cite:

- Ficklin et al. (2017) - Original KINC publication
- Shealy et al. (2019) - GPU Implementation of Pairwise GMM
- Burns et al. (2022) - Addressing Noise in Co-expression Network Construction

See :doc:`cite` for full citations.

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
