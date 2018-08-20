#!/bin/bash

set -x

KINC="build/cli/kinc"
DATA="data"
LOGS="logs"

INFILE="$1"

# import emx
EMX_FILE="$DATA/$(basename $INFILE .txt).emx"

$KINC run import-emx --input $INFILE --output $EMX_FILE --nan NA

# similarity
NP=$(nproc)
CCM_FILE="$DATA/$(basename $EMX_FILE .emx).ccm"
CMX_FILE="$DATA/$(basename $EMX_FILE .emx).cmx"
CLUSMETHOD="none"
CORRMETHOD="pearson"
MINEXPR=-inf
CRITERION="BIC"
MINCORR=0
MAXCORR=1

mpirun -np $NP $KINC run similarity --input $EMX_FILE --ccm $CCM_FILE --cmx $CMX_FILE --clusmethod $CLUSMETHOD --corrmethod $CORRMETHOD --minexpr $MINEXPR --crit $CRITERION --mincorr $MINCORR --maxcorr $MAXCORR

# threshold
LOG_FILE="$LOGS/$(basename $CMX_FILE .cmx).txt"

mkdir -p $LOGS

$KINC run rmt --input $CMX_FILE --log $LOG_FILE

# export cmx
OUTFILE="$DATA/$(basename $CMX_FILE .cmx)-cmx.txt"

$KINC run export-cmx --emx $EMX_FILE --ccm $CCM_FILE --cmx $CMX_FILE --output $OUTFILE

# extract
NET_FILE="$DATA/$(basename $EMX_FILE .emx)-net.txt"
GML_FILE="$DATA/$(basename $EMX_FILE .emx).graphml"
THRESHOLD=0

$KINC run extract --emx $EMX_FILE --ccm $CCM_FILE --cmx $CMX_FILE --output $NET_FILE --graphml $GML_FILE --mincorr $THRESHOLD
