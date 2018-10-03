#!/bin/bash

set -x

KINC="build/cli/kinc"
DATA="data"
LOGS="logs"
GPU=1

INFILE="$1"

# apply settings
if [[ $GPU == 1 ]]; then
   $KINC settings set opencl 0:0
   $KINC settings set logging off

   NP=1
else
   $KINC settings set opencl none
   $KINC settings set logging off

   NP=$(nproc)
fi

# import emx
EMX_FILE="$DATA/$(basename $INFILE .txt).emx"

$KINC run import-emx --input $INFILE --output $EMX_FILE --nan NA

# similarity
CCM_FILE="$DATA/$(basename $EMX_FILE .emx).ccm"
CMX_FILE="$DATA/$(basename $EMX_FILE .emx).cmx"
CLUSMETHOD="gmm"
CORRMETHOD="pearson"
MINEXPR="-inf"
MINCLUS=1
MAXCLUS=5
CRITERION="BIC"
PREOUT="--preout"
POSTOUT="--postout"
MINCORR=0.5
MAXCORR=1

mpirun -np $NP $KINC run similarity \
   --input $EMX_FILE \
   --ccm $CCM_FILE \
   --cmx $CMX_FILE \
   --clusmethod $CLUSMETHOD \
   --corrmethod $CORRMETHOD \
   --minexpr $MINEXPR \
   --minclus $MINCLUS --maxclus $MAXCLUS \
   --crit $CRITERION \
   $PREOUT $POSTOUT \
   --mincorr $MINCORR --maxcorr $MAXCORR

# threshold
LOG_FILE="$LOGS/$(basename $CMX_FILE .cmx).txt"

mkdir -p $LOGS

$KINC run rmt \
   --input $CMX_FILE \
   --log $LOG_FILE

# export cmx
OUTFILE="$DATA/$(basename $CMX_FILE .cmx)-cmx.txt"

$KINC run export-cmx \
   --emx $EMX_FILE \
   --ccm $CCM_FILE \
   --cmx $CMX_FILE \
   --output $OUTFILE

# extract
NET_FILE="$DATA/$(basename $EMX_FILE .emx)-net.txt"
GML_FILE="$DATA/$(basename $EMX_FILE .emx).graphml"
MINCORR=0
MAXCORR=1

$KINC run extract \
   --emx $EMX_FILE \
   --ccm $CCM_FILE \
   --cmx $CMX_FILE \
   --output $NET_FILE \
   --graphml $GML_FILE \
   --mincorr $MINCORR \
   --maxcorr $MAXCORR
