#!/bin/bash
# example: ./kinc-docker.sh cuda 4 Yeast.emx.txt

MODE="$1"
NP="$2"
INFILE="$3"

nvidia-docker run \
	--rm -it \
	-v ${PWD}:/workspace \
	systemsgenetics/kinc:latest-gpu \
	bash -c "./kinc.sh ${MODE} ${NP} ${INFILE}; chown $(stat -c '%u:%g' .) *"
