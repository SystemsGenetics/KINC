[![DOI](https://zenodo.org/badge/71836133.svg)](https://zenodo.org/badge/latestdoi/71836133)

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

To build & install KINC:
```
cd build
qmake ../src/KINC.pro
make qmake_all
make
make qmake_all
make install
```

## Using the KINC GUI or Console

ACE provides two different libraries for GUI and console applications. The `kinc` executable is the console or command line version and the `qkinc` executable is the GUI version.

# Usage

To build a GCN involves several steps:

1. Import expression matrix
2. Compute cluster composition matrix
3. Compute correlation matrix
4. Compute thresholded correlation matrix

# Troubleshooting
## An error occurred in MPI_Init
KINC requires MPI as a dependency, but on most systems you can execute the command-line KINC as a stand-alone tool without using 'mpirun'.  This is because KINC checks during runtime if MPI is appropriate for execution. However, on a SLURM cluster where MPI jobs must be run using the srun command and where PMI2 is compiled into MPI, then KINC cannot be executed stand-alone.  It must be executed using srun with the --mpi argument set to pmi2.  For example:

```
srun --mpi=pmi2 kinc run import_emx --input Yeast-ematrix.txt --output Yeast.emx --nan NA
```

