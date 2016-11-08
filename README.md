RMTGeneNet
==========

This package provides two binaries. One that generates pearson correlations for
gene expression data and another that employs Random Matrix Theory (RMT) to 
determine thresholding of a gene co-expression network. A perl script can be
used to extract the gene co-expression network using the RMT threshold from
the correlation binary file.

This package was converted from Java source code to C code to improve 
performance.  The original Java code can be found here:
http://bci.clemson.edu/software/rmt.

RMTGeneNet is offered under a GPL v2.0 license agreement.  


INSTALLATION
------------
For installation instructions see the INSTALL.txt file.


BUILD A NETWORK
---------------
To construct a gene co-expression network using RMTGeneNet an expression
matrix must first be constructed.  The expression matrix should be present in
a tab-delimeted file where the columns represent the experimental samples and
the rows represent the measured unit. In the case of a set of microarray
experiments the columns would represent each microarray and the rows would
represent the probesets of the array.  The first column should contain the
names of measured units (e.g. probesets) but there should be no column headers
in the file.  The cells in the matrix should contain the measured expression 
levels. The file should be named with a .txt extension.

Step 1: Calculate Pearson Corrleations
Use the 'ccm' binary to construct pearson correlations.  The binary can be
executed in the following way:

   ccm <ematrix> <rows> <cols> [<omit_na> <na_val> <hist> <perf> <headers>]

Where:
    <ematrix>: the file name that contains the expression matrix. The rows must 
               be genes or probesets and columns are samples.
    <rows>:    the number of lines in the input file minus the header column if 
               it exists.
    <cols>:    the number of columns in the input file minus the first column 
               that contains gene names.
    <omit_na>: set to 1 to ignore missing values. Defaults to 0.
    <na_val>:  a string representing the missing values in the input file 
               (e.g. NA or 0.000).
    <min_obs>: the minimum number of observations (after missing values removed) 
               that must be present to perform correlation. Default is 30.
    <func>:    a transformation function to apply to elements of the ematrix. 
               Values include: log2 or none. Default is none.
    <hist>:    set to 1 to enable creation of correlation historgram. Defaults 
               to 0.
    <perf>:    set to 1 to enable performance monitoring. Defaults to 0.
    <headers>: set to 1 if the first line contains headers. Defaults to 0.
    

The 'ccm' binary will create several binary files inside of a 'Pearson'
directory.  These binary files contain the pair-wise Pearson correlation
values for every unit.  Correlation value is set to NaN if there weren't enough 
observations to perform the calculation.


Step 2:  Calculate a Network Threshold
After the Pearson correlations have been calculated, the 'rmm' employs Random
Matrix Theory to identify a suitable threshold where the nearest neighbor
spacing distribution of the eigenvalues.  See the following reference for
further details:

Luo F, Yang Y, Zhong J, Gao H, Khan L, Thompson DK, Zhou J (2007) Constructing
   gene co-expression networks and predicting functions of unknown genes by 
   random matrix theory.  BMC Bioinformatics 8: 299

To find the network threshold execute the 'rmm' with the following arguments:

    '-i': The input file name. Same as used in previous step.
        Must be the same as the name used in the matrix binary
        files.  This name will also be used to create files 
        associated with this run of the program.
        
Optional:
    
    '-b': The initial threshold(+1*step) value that will be used.
        [b=Begin value] Default: 0.9200
    '-s': The threshold step size used each iteration. Default: 0.001
    '-c': The chi-square test value that the loop will stop on.
        Default: 200
    '-v': Set the performance collection. Has two values possible values,
        ON/OFF . [v=Verbose] Default: ON
        
Examples:

    <executable> -i <input.file.name> 
    <exec> -i <input.file.name> -s <0.0001> -v ON
    
Argument order is not important, but spaces are required between the flag and  
value. If the underlying data set is highly correlated, the 'rmm' binary may
not be able to find a threshold.

Step 3: Generate the Network
Once a threshold has been determined, the Perl script named
'parse_pearson_bin.pl' and be used to generate the network file.  Execute the
script in the following way:

   parse_pearson_bin.pl -b Pearson -t <threshold> -p <probesets> -o <prefix>

Where

   <threshold> is the numeric threshold value provided by the 'rmm' binary
   <probestes> is a text file containing the measurable units (e.g.
               probesets). They must appear in the same order, one unit per 
               line, as they appear in the rows of the expression matrix.
   <prefix>    a prefix to attach to the output network name.

The '-b Pearson' argument specifies the Pearson directory created in step1
above where the Pearson correlation binary files are located.  

The 'parse_pearson_bin.pl' script has other functionality as well. For further
details, execute the script with a '--help' argument.




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


