# Knowledge Independent Network Construction (KINC)

KINC is a Qt/ACE application that produces a gene co-expression network (GCN) from a gene expression matrix (GEM). KINC implements several algorithms that are required for this process:

Correlation
- Pearson
- Spearman

Clustering
- Gaussian mixture models

Thresholding
- Power-law
- Random matrix theory

KINC is built with [ACE](https://github.com/SystemsGenetics/ACE), a framework which provides mechanisms for large-scale heterogeneous computing and data management. As such, KINC can be run in a variety of compute configurations, including single-core / single-GPU and multi-core / multi-GPU, and KINC uses its own binary file formats to represent the data objects that it produces. Each of these binary formats can be exported to a plain-text format for use in other applications.

## Installation

Refer to the files under `docs` for installation instructions. KINC is currently supported on most flavors of Linux.

## Usage

KINC provides two executables: `kinc`, the command-line version, and `qkinc`, the GUI version. The command-line version can use MPI while the GUI version can display data object files that are produced by KINC. KINC produces a gene-coexpression network in several steps:
1. `import-emx`: Import expression matrix text file into binary format
2. `similarity`: Compute a cluster matrix and correlation matrix from expression matrix
3. `threshold`: Determine an appropriate correlation threshold for correlation matrix
4. `extract`: Extract an edge list from a correlation matrix given a threshold

Below is an example usage of `kinc` on the Yeast dataset:
```
# import expression matrix into binary format
kinc run import-emx --input Yeast-GEM.txt --output Yeast.emx --nan NA

# compute similarity matrix (with GMM clustering)
mpirun -np 8 kinc run similarity --input Yeast.emx --ccm Yeast.ccm --cmx Yeast.cmx --clusmethod gmm --corrmethod spearman --minclus 1 --maxclus 5

# determine correlation threshold
kinc run rmt --input Yeast.cmx --log Yeast.log

# read threshold from log file
THRESHOLD=$(tail -n 1 Yeast.log)

# extract network file from thresholded similarity matrix
kinc run extract --emx Yeast.emx --ccm Yeast.ccm --cmx Yeast.cmx --output Yeast-net.txt --mincorr $THRESHOLD
```

A more thorough example usage is provided in `scripts/run-all.sh`.

### Running KINC on SLURM

Although KINC is an MPI application, generally you can run `kinc` as a stand-alone application without `mpirun` and achieve normal serial behavior. However, on a SLURM cluster where MPI jobs must be run with the `srun` command and where PMI2 is compiled into MPI, `kinc` cannot be executed stand-alone. It must be executed using `srun` with the additional argument `--mpi=pmi2`. For example:
```
srun --mpi=pmi2 kinc run import_emx --input Yeast-ematrix.txt --output Yeast.emx --nan NA
```
