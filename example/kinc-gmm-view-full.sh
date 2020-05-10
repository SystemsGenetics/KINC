#!/bin/bash
# View the entire network

kinc-3d-viewer.py \
    --net "PRJNA301554.slim.GEM.log2.paf-th0.00-p1e-3-rsqr0.30-filtered-th_ranked.csGCN.txt" \
    --emx "PRJNA301554.slim.GEM.log2.txt" \
    --amx "PRJNA301554.slim.annotations.txt"
