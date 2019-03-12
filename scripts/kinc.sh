#!/bin/bash

# parse command-line arguments
if [[ $# != 2 ]]; then
	echo "usage: $0 <mode> <infile>"
	exit -1
fi

# define execution mode
MODE="$1"

# define analytic flags
DO_IMPORT_EMX=1
DO_SIMILARITY=1
DO_EXPORT_CMX=0
DO_THRESHOLD=0
DO_EXTRACT=0

# define input/output files
INFILE="$2"
DATA="data"
EMX_FILE="$DATA/$(basename $INFILE .txt).emx"
CCM_FILE="$DATA/$(basename $EMX_FILE .emx).ccm"
CMX_FILE="$DATA/$(basename $EMX_FILE .emx).cmx"
LOGS="logs"
RMT_FILE="$LOGS/$(basename $CMX_FILE .cmx).txt"

# apply settings
if [[ $MODE = "cuda" ]]; then
	kinc settings set cuda 0
	kinc settings set opencl none
	kinc settings set threads 2
	kinc settings set logging off

	NP=1
elif [[ $MODE = "opencl" ]]; then
	kinc settings set cuda none
	kinc settings set opencl 0:0
	kinc settings set threads 2
	kinc settings set logging off

	NP=1
elif [[ $MODE == "serial" ]]; then
	kinc settings set cuda none
	kinc settings set opencl none
	kinc settings set logging off

	NP=$(nproc)
else
	echo "error: invalid execution mode \'$MODE\'"
	exit -1
fi

# import emx
if [[ $DO_IMPORT_EMX = 1 ]]; then
	time kinc run import-emx \
		--input $INFILE \
		--output $EMX_FILE \
		--nan NA
fi

# similarity
if [[ $DO_SIMILARITY = 1 ]]; then
	CLUSMETHOD="gmm"
	CORRMETHOD="spearman"
	MINEXPR="-inf"
	MINCLUS=1
	MAXCLUS=5
	CRITERION="BIC"
	PREOUT="--preout"
	POSTOUT="--postout"
	MINCORR=0.5
	MAXCORR=1
	LSIZE=32

	time mpirun -np $NP kinc run similarity \
		--input $EMX_FILE \
		--ccm $CCM_FILE \
		--cmx $CMX_FILE \
		--clusmethod $CLUSMETHOD \
		--corrmethod $CORRMETHOD \
		--minexpr $MINEXPR \
		--minclus $MINCLUS \
		--maxclus $MAXCLUS \
		--crit $CRITERION \
		$PREOUT $POSTOUT \
		--mincorr $MINCORR \
		--maxcorr $MAXCORR \
		--lsize $LSIZE
fi

# export cmx
if [[ $DO_EXPORT_CMX = 1 ]]; then
	OUTFILE="$DATA/$(basename $CMX_FILE .cmx)-cmx.txt"

	time kinc run export-cmx \
		--emx $EMX_FILE \
		--ccm $CCM_FILE \
		--cmx $CMX_FILE \
		--output $OUTFILE
fi

# threshold
if [[ $DO_THRESHOLD = 1 ]]; then
	mkdir -p $LOGS

	time kinc run rmt \
	   --input $CMX_FILE \
	   --log $RMT_FILE
fi

# extract
if [[ $DO_EXTRACT = 1 ]]; then
	NET_FILE="$DATA/$(basename $EMX_FILE .emx)-net.txt"
	MINCORR=0
	MAXCORR=1

	time kinc run extract \
		--emx $EMX_FILE \
		--ccm $CCM_FILE \
		--cmx $CMX_FILE \
		--output $NET_FILE \
		--mincorr $MINCORR \
		--maxcorr $MAXCORR
fi
