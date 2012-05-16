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
First create a 'Pearson' directory:

   mkdir Pearson

Next, use the 'ccm' binary to construct pearson correlations.  The binary can be
executed in the following way:

   ccm <input file> <rows> <cols> 

Where:
   <input file> is the name of the expression matrix file.  Do not include the
                .txt extension, just the name.
   <rows>       the number of measurable units (or rows) in the matrix
   <cols>       the number of samples in the matrix. This should be one minus
                the actual number of columns because the first column is the
                probeset name.

The 'ccm' binary will create several binary files inside of the 'Pearson'
directory.  These binary files contain the pair-wise Pearson correlation
values for every unit.


Step 2:  Calculate a Network Threshold
After the Pearson correlations have been calculated, the 'rmm' employs Random
Matrix Theory to identify a suitable threshold where the nearest neighbor
spacing distribution of the eigenvalues.  See the following reference for
further details:

Luo F, Yang Y, Zhong J, Gao H, Khan L, Thompson DK, Zhou J (2007) Constructing
   gene co-expression networks and predicting functions of unknown genes by 
   random matrix theory.  BMC Bioinformatics 8: 299

To find the network threshold execute the 'rmm' binary in the following way:

   rmm -i <input file> 

Where:

   <input file> is the name of the expression matrix file.  Do not include the
                .txt extension, just the name.

By default the 'rmm' binary will being selection of the threshold at a
correlation value of 0.92. If this is not the correct threshold, it will
repeatedly decrease the correlation value in increments of 0.001 until it
finds the threshold.  The starting threshold and step increment can be
specified at run time.  For further instructions, execute the 'rmm' binary
with no arguments.

Note:  If the underlying data set is highly correlated, the 'rmm' binary may
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




