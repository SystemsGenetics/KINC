#!/bin/bash
# View unique Duration edges (i.e. not covariant with Treatment/Genotype/Subspecies)

kinc-3d-viewer.py \
    --net "PRJNA301554.slim.GEM.log2.paf-th0.00-p1e-3-rsqr0.30-filtered-th_ranked.Duration-unique_class.csGCN.txt" \
    --emx "../data/PRJNA301554.slim.GEM.log2.txt" \
    --amx "../data/PRJNA301554.slim.annotations.txt"
