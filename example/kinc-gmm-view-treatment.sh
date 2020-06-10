#!/bin/bash
# View unique Treatment edges (i.e. not covariant with Duration/Genotype/Subspecies)

kinc-3d-viewer.py \
    --net "PRJNA301554.slim.GEM.log2.paf-th0.00-p1e-3-rsqr0.30-filtered-th_ranked.Treatment-unique_class.csGCN.txt" \
    --emx "../data/PRJNA301554.slim.GEM.log2.txt" \
    --amx "../data/PRJNA301554.slim.annotations.txt"
