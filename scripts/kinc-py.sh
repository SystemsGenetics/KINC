#!/bin/bash

# parse command-line arguments
if [[ $# != 1 ]]; then
	echo "usage: $0 <infile>"
	exit -1
fi

# define analytic flags
DO_SIMILARITY=1
DO_THRESHOLD=1
DO_EXTRACT=1

# define input/output files
INFILE="$1"
DIRNAME="$(dirname ${INFILE})"
BASENAME="$(basename ${INFILE} .txt)"
EMX_FILE="${INFILE}"
CMX_FILE="${DIRNAME}/${BASENAME}.cmx-py.txt"

# similarity
if [[ ${DO_SIMILARITY} = 1 ]]; then
	CLUSMETHOD="gmm"
	CORRMETHOD="spearman"
	MINEXPR="-inf"
	MINCLUS=1
	MAXCLUS=5
	CRITERION="bic"
	PREOUT="--preout"
	POSTOUT="--postout"
	MINCORR=0.5
	MAXCORR=1

	env time -f "%e" python scripts/kinc-similarity.py \
		--input ${EMX_FILE} \
		--output ${CMX_FILE} \
		--clusmethod ${CLUSMETHOD} \
		--corrmethod ${CORRMETHOD} \
		--minexpr=${MINEXPR} \
		--minclus ${MINCLUS} \
		--maxclus ${MAXCLUS} \
		--criterion ${CRITERION} \
		${PREOUT} \
		${POSTOUT} \
		--mincorr ${MINCORR} \
		--maxcorr ${MAXCORR}
fi

# threshold
if [[ ${DO_THRESHOLD} = 1 ]]; then
	NUM_GENES=$(expr $(cat ${EMX_FILE} | wc -l) - 1)
	METHOD="rmt"
	TSTART=0.99
	TSTEP=0.001
	TSTOP=0.50

	env time -f "%e" python scripts/kinc-threshold.py \
		--input ${CMX_FILE} \
		--n-genes ${NUM_GENES} \
		--method ${METHOD} \
		--tstart ${TSTART} \
		--tstep ${TSTEP} \
		--tstop ${TSTOP}
fi

# extract
if [[ ${DO_EXTRACT} = 1 ]]; then
	MINCORR=0
	MAXCORR=1
	NET_FILE="${DIRNAME}/${BASENAME}.th000.coexpnet-py.txt"

	env time -f "%e" python scripts/kinc-extract.py \
		--emx ${EMX_FILE} \
		--cmx ${CMX_FILE} \
		--output ${NET_FILE} \
		--mincorr ${MINCORR} \
		--maxcorr ${MAXCORR}
fi
