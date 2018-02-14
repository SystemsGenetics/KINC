
![alt tag](https://raw.githubusercontent.com/SystemsGenetics/KINC/version1/KINClogo.png)

# Knowledge Independent Network Construction

KINC is an application that produces a gene co-expression network (GCN) from a gene expression matrix (GEM).

## Installation

KINC depends on the following software:

- GNU Scientific Library (http://www.gnu.org/software/gsl/)
- mixmod (http://www.mixmod.org/)
- BLAS (http://www.netlib.org/blas/)
- LAPACK (http://www.netlib.org/lapack/)

To build KINC:
```
make
make install [INSTALL_PREFIX=...]
```

## Usage

KINC produces a GCN in three steps: `similarity`, `threshold`, and `extract`. Run `kinc help` to view usage for `kinc` and its sub-commands.

Below is an example usage of KINC on the Yeast dataset:
```
# setup input data
GEMFILE="Yeast-GEM.txt"
ROWS=$(wc -l < $GEMFILE)
COLS=$(awk -F ' ' '{print NF ; exit}' $GEMFILE)

# compute similarity matrix (with GMM clustering)
./kinc similarity --ematrix $GEMFILE --rows $ROWS --cols $COLS --headers --clustering mixmod --criterion BIC --method sc --min_obs 30 --omit_na --na_val NA

# determine correlation threshold
./kinc threshold --ematrix $GEMFILE --rows $ROWS --cols $COLS --headers --clustering mixmod --method sc --th_method sc --max_modes 5 --omit_na --na_val NA

# extract network file from thresholded similarity matrix
THRESHOLD=$(cat Yeast-GEM.sc.mcs30.md5.mmINF.th.txt)

./kinc extract --ematrix $GEMFILE --rows $ROWS --cols $COLS --headers --clustering mixmod --method sc --th_method sc --max_modes 5 --omit_na --na_val NA --th $THRESHOLD
```

### Input

KINC takes a gene expression matrix as input. The expression matrix should be a tab-delimited file where the columns represent the experimental samples and the rows represent the measured unit. In the case of a set of microarray experiments the columns would represent each microarray and the rows would represent the probesets of the array.  For RNA-Seq, the columns represent individual samples and the columns represent genes or transcripts. The first column should contain the names of measured units (e.g. probesets or genes) and the first row may contain sample names. The cells in the matrix should contain the measured expression levels.

### Similarity

Compute a similarity matrix of all genes using Pearson or Spearman correlation. Gaussian mixture model (GMM) clustering can be optionally employed to identify different modes of interaction in each gene pair. However, note that clustering significantly increases the runtime and the size of the output data.

### Threshold

Determine the appropriate threshold to apply to the similarity matrix in order to create an adjacency matrix. KINC finds this threshold by starting at the maximum and decrementing until the appropriate threshold is found using random matrix theory (RMT):

1. Compute thresholded similarity matrix (cut matrix)
2. Compute eigenvalues of the cut matrix
3. Compute the nearest neighbor spacing distribution of the eigenvalues
4. Stop when the distribution transitions from Poisson to Gaussian

For further information on this technique refer to the following paper:

```
Luo F, Yang Y, Zhong J, Gao H, Khan L, Thompson DK, Zhou J (2007) Constructing
   gene co-expression networks and predicting functions of unknown genes by
   random matrix theory.  BMC Bioinformatics 8: 299
```

### Extract

Given a correlation threshold, extract the relevant correlations from the similarity matrix. If clustering was used during similarity, the output will also contain cluster information.
