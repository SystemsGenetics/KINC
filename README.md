First, open a new KINC Expression Matrix file 

open Hsapiens-9606-201603-2016-RNASeq-Quantile-CancerGenomeAtlas-v1.emx:emx --select

Second, import the file the KINC file the expression data from the flat text file. This file has headers and empty values are indicated with 'NA'

load Hsapiens-9606-201603-2016-RNASeq-Quantile-CancerGenomeAtlas-v1 --nosample=NA

Next see what GPU devices are available
cl list

Set the CL device:

cl set 0:0

Now run spearman
spearman --in=Hsapiens-9606-201603-2016-RNASeq-Quantile-CancerGenomeAtlas-v1.emx --out=Hsapiens-9606-201603-2016-RNASeq-Quantile-CancerGenomeAtlas-v1.cmx:cmx
