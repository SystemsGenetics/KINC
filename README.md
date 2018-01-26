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
sudo apt install libgsl2 ocl-icd-opencl-dev libopenmpi-dev
```

To build KINC:
```
cd build
qmake ../src
make
cd ..
```

## Using the KINC GUI or Console

ACE provides two different libraries for GUI and console applications. The `GUI` variable in `KINC.pro` controls which library KINC uses. By default, KINC is built as a GUI application. When run as a GUI application, KINC prints the console equivalent of every analytic as it is run.

# Usage

To build a GCN involves several steps:

1. Import expression matrix
2. Compute cluster composition matrix
3. Compute correlation matrix
4. Compute thresholded correlation matrix
