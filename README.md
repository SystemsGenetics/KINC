
![alt tag](https://raw.githubusercontent.com/SystemsGenetics/KINC/version1/KINClogo.png)

Knowledge Indpendent Network Construction
==========

This package provides two binaries. One that generates pearson correlations for
gene expression data and another that employs Random Matrix Theory (RMT) to 
determine thresholding of a gene co-expression network. A perl script can be
used to extract the gene co-expression network using the RMT threshold from
the correlation binary file.

This package was converted from Java source code to C code to improve 
performance.  The original Java code can be found here:
http://bci.clemson.edu/software/rmt.

KINC is offered under a GPL v2.0 license agreement.  


INSTALLATION
------------
For installation instructions see the INSTALL.txt file.


BUILD A NETWORK
---------------
To construct a gene co-expression network using KINC, an expression
matrix must first be constructed.  The expression matrix should be present in
a tab-delimeted file where the columns represent the experimental samples and
the rows represent the measured unit. In the case of a set of microarray
experiments the columns would represent each microarray and the rows would
represent the probesets of the array.  For RNA-Seq, the columns represent individual
samples and the columns represent genes or transcripts.  The first column 
should contain the names of measured units (e.g. probesets or genes) but there 
should be no column headers in the file.  The cells in the matrix should contain 
the measured expression levels. 

Step 1: Calculate Similarity Matrix (Correlation)

TODO: add instructions

Step 2:  Determine a Threshold
After the correlations have been calculated a threshold should be determined
to conver the simialarity matrix into an adjacency matrix.  Actually, KINC
does not create an adjacency matrix, but the threshold to do so must be
found.  KINC currently employes Random Matrix Theory (RMT) to identify a suitable 
threshold where the nearest neighbor spacing distribution of the eigenvalues ceases
to appear Poisson.  See the following reference for further details:

Luo F, Yang Y, Zhong J, Gao H, Khan L, Thompson DK, Zhou J (2007) Constructing
   gene co-expression networks and predicting functions of unknown genes by 
   random matrix theory.  BMC Bioinformatics 8: 299

To find the network threshold execute the 'rmm' with the following arguments:

TODO: add instructions

Step 3: Extract the Network
Once a threshold has been determined we can generate flat text files that 
represent the network.  


TODO: add instructions



