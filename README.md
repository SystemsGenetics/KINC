# Knowledge Independent Network Construction (KINC)

KINC is a Qt/ACE application that produces a gene co-expression network (GCN) from a gene expression matrix (GEM). KINC implements several algorithms that are required for this process:

Correlation
- Pearson
- Spearman

Clustering
- K-means
- Gaussian mixture models

Thresholding
- Random matrix theory

# Installation

This software uses GSL, OpenCL, and [ACE](https://github.com/SystemsGenetics/ACE). For instructions on installing ACE, see the project repository. For all other dependencies, consult your package manager. For example, to install dependencies on Ubuntu:
```
sudo apt install libgsl2 ocl-icd-opencl-dev
```

To build KINC:
```
cd build
qmake ../src
make
cd ..
```

# Usage

To build a GCN involves several steps:

1. Import expression matrix
2. Compute cluster composition matrix
3. Compute correlation matrix
4. Compute thresholded correlation matrix

KINC provides a GUI to perform each step one at a time, as well as a console to show the equivalent command-line operation. That is, every KINC operation can be done through the GUI or the command-line.
