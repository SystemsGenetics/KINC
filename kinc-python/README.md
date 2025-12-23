# KINC-Python

Python implementation of Knowledge Independent Network Construction (KINC) for constructing gene co-expression networks from gene expression matrices.

## Overview

KINC-Python is a complete rewrite of KINC in Python, providing:
- GPU-accelerated pairwise correlation computation (Pearson, Spearman)
- Gaussian Mixture Model (GMM) based clustering for condition-specific networks
- MPI support for distributed computing on HPC clusters
- Multi-core CPU parallelization for workstations
- HDF5-based data formats for efficiency and interoperability

## Features

- **Core Analytics**:
  - `import-emx`: Import gene expression matrices from text files
  - `export-emx`: Export expression matrices to text
  - `similarity`: Compute pairwise correlations with optional GMM clustering
  - `corrpower`: Filter edges by statistical power
  - `cond-test`: Test for condition-specific correlations
  - `rmt`: Random Matrix Theory thresholding
  - `powerlaw`: Power-law thresholding
  - `extract`: Extract network edges to text/GraphML

- **Parallelization**:
  - MPI distribution across compute nodes (mpi4py)
  - GPU acceleration with CuPy/RAPIDS
  - Multi-core CPU parallelization with Dask

- **Data Formats**:
  - HDF5 (.h5) for efficient storage and fast I/O
  - Compatible with standard Python scientific tools
  - Self-documenting with embedded metadata

## Installation

```bash
# Basic installation
pip install -e .

# With GPU support
pip install -e ".[gpu]"

# With MPI support
pip install -e ".[mpi]"

# Full installation (all features)
pip install -e ".[all]"
```

## Requirements

- Python 3.8+
- NumPy, SciPy, pandas
- h5py (HDF5 support)
- scikit-learn (GMM clustering)
- Optional: cupy, mpi4py, dask

## Quick Start

```bash
# Import gene expression matrix
kinc import-emx --input data.txt --output data.emx.h5

# Compute similarity with GMM clustering
kinc similarity \
  --input data.emx.h5 \
  --output data.cmx.h5 \
  --ccm data.ccm.h5 \
  --clusmethod gmm \
  --corrmethod spearman \
  --minclus 1 \
  --maxclus 5

# Extract network
kinc extract \
  --input data.cmx.h5 \
  --output network.txt \
  --mincorr 0.85
```

## Architecture

```
kinc-python/
├── kinc/
│   ├── core/          # Core algorithms (correlation, GMM, etc.)
│   ├── data/          # Data structures (ExpressionMatrix, etc.)
│   ├── analytics/     # Analytics (similarity, rmt, etc.)
│   ├── parallel/      # Parallelization backends
│   └── cli.py         # Command-line interface
├── tests/             # Unit tests
├── examples/          # Example workflows
└── docs/              # Documentation
```

## Citation

If you use KINC-Python in your research, please cite:

- Ficklin et al. (2017) - Original KINC publication
- Shealy et al. (2019) - GPU Implementation of Pairwise GMM
- Burns et al. (2022) - Addressing Noise in Co-expression Network Construction

## License

This project maintains compatibility with the original KINC license.

## Development Status

This is a faithful Python reimplementation of KINC v3.4.2, designed to be statistically equivalent to the C++ version while leveraging modern Python scientific libraries.
