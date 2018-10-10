#!/bin/bash

# parse command-line arguments
if [[ $# != 1 ]]; then
	echo "usage: $0 <infile>"
	exit -1
fi

DATA="data"
EMX_FILE="$1"

# similarity
CMX_FILE="$DATA/$(basename $EMX_FILE .txt)-cmx-py.txt"
CLUSMETHOD="gmm"
CORRMETHOD="pearson"
MINEXPR="-inf"
MINCLUS=1
MAXCLUS=5
CRITERION="bic"
PREOUT="--preout"
POSTOUT="--postout"
MINCORR=0.5
MAXCORR=1

python scripts/similarity.py \
   -i $EMX_FILE \
   -o $CMX_FILE \
   --clusmethod $CLUSMETHOD \
   --corrmethod $CORRMETHOD \
   --minexpr=$MINEXPR \
   --minclus $MINCLUS --maxclus $MAXCLUS \
   --crit $CRITERION \
   $PREOUT $POSTOUT \
   --mincorr $MINCORR --maxcorr $MAXCORR

# threshold
NUM_GENES=$(expr $(cat $EMX_FILE | wc -l) - 1)
METHOD="rmt"
TSTART=0.99
TSTEP=0.001
TSTOP=0.50

python scripts/threshold.py \
   -i $CMX_FILE \
   --genes $NUM_GENES \
   --method $METHOD \
   --tstart $TSTART \
   --tstep $TSTEP \
   --tstop $TSTOP

# extract
NET_FILE="$DATA/$(basename $EMX_FILE .txt)-net-py.txt"
MINCORR=0
MAXCORR=1

python scripts/extract.py \
   --emx $EMX_FILE \
   --cmx $CMX_FILE \
   --output $NET_FILE \
   --mincorr $MINCORR \
   --maxcorr $MAXCORR
