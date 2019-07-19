#!/bin/bash

# parse command-line arguments
if [[ $# != 3 ]]; then
	echo "usage: $0 <mode> <np> <infile>"
	exit -1
fi

# define execution mode
MODE="$1"
NP="$2"

# define analytic flags
DO_IMPORT_EMX=1
DO_SIMILARITY=1
DO_CORRPOWER=0
DO_EXPORT_CMX=0
DO_THRESHOLD=1
DO_EXTRACT=1

# define input/output files
INFILE="$3"
DIRNAME="$(dirname ${INFILE})"
BASENAME="$(basename ${INFILE} .txt)"
EMX_FILE="${DIRNAME}/${BASENAME}.emx"
CCM_FILE="${DIRNAME}/${BASENAME}.ccm"
CMX_FILE="${DIRNAME}/${BASENAME}.cmx"
CMX_TXT_FILE="${DIRNAME}/${BASENAME}.cmx.txt"
RMT_FILE="${DIRNAME}/${BASENAME}.rmt.txt"

# apply settings
kinc settings set threads 2
kinc settings set logging off

if [[ ${MODE} = "cuda" ]]; then
	kinc settings set cuda 0
	kinc settings set opencl none
elif [[ ${MODE} = "opencl" ]]; then
	kinc settings set cuda none
	kinc settings set opencl 0:0
elif [[ ${MODE} == "serial" ]]; then
	kinc settings set cuda none
	kinc settings set opencl none
else
	echo "error: invalid execution mode \'${MODE}\'"
	exit -1
fi

# import emx
if [[ ${DO_IMPORT_EMX} = 1 ]]; then
	env time -f "%e" kinc run import-emx \
		--input ${INFILE} \
		--output ${EMX_FILE} \
		--nan NA
fi

# similarity
if [[ ${DO_SIMILARITY} = 1 ]]; then
	CLUSMETHOD="gmm"
	CORRMETHOD="spearman"
	MINEXPR="-inf"
	MINCLUS=1
	MAXCLUS=5
	CRITERION="ICL"
	PREOUT="true"
	POSTOUT="true"
	MINCORR=0.5
	MAXCORR=1
	LSIZE=32

	env time -f "%e" mpirun -np ${NP} kinc run similarity \
		--input ${EMX_FILE} \
		--ccm ${CCM_FILE} \
		--cmx ${CMX_FILE} \
		--clusmethod ${CLUSMETHOD} \
		--corrmethod ${CORRMETHOD} \
		--minexpr ${MINEXPR} \
		--minclus ${MINCLUS} \
		--maxclus ${MAXCLUS} \
		--crit ${CRITERION} \
		--preout ${PREOUT} \
		--postout ${POSTOUT} \
		--mincorr ${MINCORR} \
		--maxcorr ${MAXCORR} \
		--lsize ${LSIZE}
fi

# correlation power analysis
if [[ ${DO_CORRPOWER} = 1 ]]; then
	CCM_OUT_FILE="${DIRNAME}/${BASENAME}.corrpower.ccm"
	CMX_OUT_FILE="${DIRNAME}/${BASENAME}.corrpower.cmx"
	ALPHA=0.001
	POWER=0.8

	env time -f "%e" mpirun -np ${NP} kinc run corrpower \
		--ccm-in ${CCM_FILE} \
		--cmx-in ${CMX_FILE} \
		--ccm-out ${CCM_OUT_FILE} \
		--cmx-out ${CMX_OUT_FILE} \
		--alpha ${ALPHA} \
		--power ${POWER}
fi

# export cmx
if [[ ${DO_EXPORT_CMX} = 1 ]]; then
	env time -f "%e" kinc run export-cmx \
		--emx ${EMX_FILE} \
		--ccm ${CCM_FILE} \
		--cmx ${CMX_FILE} \
		--output ${CMX_TXT_FILE}
fi

# threshold
if [[ ${DO_THRESHOLD} = 1 ]]; then
	env time -f "%e" kinc run rmt \
		--input ${CMX_FILE} \
		--log ${RMT_FILE} \
		--threads ${NP}
fi

# extract
if [[ ${DO_EXTRACT} = 1 ]]; then
	MINCORR=0
	MAXCORR=1
	NET_FILE="${DIRNAME}/${BASENAME}.th000.coexpnet.txt"

	env time -f "%e" kinc run extract \
		--emx ${EMX_FILE} \
		--ccm ${CCM_FILE} \
		--cmx ${CMX_FILE} \
		--output ${NET_FILE} \
		--mincorr ${MINCORR} \
		--maxcorr ${MAXCORR}
fi
