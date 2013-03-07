---------------------------------------------------------
Example: Construct an S. cerevisiae (yeast) Global Network
---------------------------------------------------------
This examples provides datasets and instructions for constructing
a global co-expression networking using RMTGeneNet from microarray
samples using the Affymetrix Yeast Genome 2.0 Array. These samples
are publicly available in NCBI GEO. A list of these samples can be
found the in the file named 'samples.txt' in this directory. In
total, expression levels from 1535 samples are included in this dataset.

For this example, the file 'yeast-s_cerevisiae1.global.RMA.nc-no-na-nh.txt'
contains the expression matrix for all 1535 samples.  The columns of the
matrix correspond to the samples and the rows correspond to the probesets
on the array.  The column headers have been removed from this file, but
the exact order of samples for the columns is found in the 'samples.txt' 
file. The expression levels in this file have already been
normalized.  Control probes and ambiguous probes have been removed as well
as outlier samples. In total, 10359 probesets remain in the file.  For 
detailed instructions about the quality control steps used, please see the 
publication:

Gibson SM, Ficklin SP, Isaacson S, Luo F, Feltus FA, et al. (2013)
Massive-Scale Gene Co-Expression Network Construction and Robustness Testing
Using Random Matrix Theory. PLoS ONE 8(2): e55871.


Step 1: Construct correlation matrix
------------------------------------
The first step in construction of the yeast global network is to construct
the correlation matrix.  The following commands can be executed to construct
this matrix:

  mkdir Pearson
  ccm yeast-s_cerevisiae1.global.RMA.nc-no-na-nh 10359 1535

The number of probesets and samples in the file are passed to 'ccm'.


Step 2: Use RMT to determine an appropriate threshold
-----------------------------------------------------
The second step is to use Random Matrix Theory (RMT) to identify an
appropriate threshold for the global network. This is performed 
using the 'rmm' executable:

rmm -i yeast-s_cerevisiae1.global.RMA.nc-no-na-nh -b 0.930000

To speed up execution time we will start examining the threshold at a level
of 0.930000. Normally, you would not know where to begin looking for the
thrsehold but because this network has been build before we know it is
safe to start at this threshold.


Step 3: Generate additional network files
-----------------------------------------
The threshold returned from Step 2 was 0.831100. We can now use that
threshold to generate the final network file

perl parse_pearson_bin.pl -b Pearson -t 0.831100 -p probest_order.txt -o yeast-s_cerevisiae1.global


