#!/bin/bash
# Perform a traditional run of KINC

# Step 1: Import the similarity matrix
echo "Importing the gene expression matrix (GEM) for the slimmed experiment PRJNA301554"

kinc run import-emx \
   --input "PRJNA301554.slim.GEM.log2.txt" \
   --output "PRJNA301554.slim.GEM.log2.emx" \
   --nan "NA" \
   --samples 0

# Step 2: Perform the similarity calculations
echo "Performing similarity calculations using traditional correlation only."

kinc run similarity \
  --input "PRJNA301554.slim.GEM.log2.emx" \
  --ccm "PRJNA301554.slim.GEM.log2.traditional.ccm" \
  --cmx "PRJNA301554.slim.GEM.log2.traditional.cmx" \
  --clusmethod "none" \
  --corrmethod "spearman" \
  --minsamp 30 \
  --minexpr -inf \
  --mincorr 0.5 \
  --maxcorr 1

# Step 3:  Use RMT to find the threshold
echo "Performing threshold calculation using Random Matrix Theory (RMT)."

kinc run rmt \
  --input "PRJNA301554.slim.GEM.log2.traditional.cmx" \
  --log "PRJNA301554.slim.GEM.log2.traditional.rmt.log" \
  --tstart "0.95" \
  --tstep "0.001" \
  --tstop "0.5" \
  --threads 1 \
  --epsilon 1e-6 \
  --mineigens 50 \
  --spline true \
  --minpace 10 \
  --maxpace 40 \
  --bins 60

# Step 5: Extract the network using the RMT threshold.
echo "Extracting the RMT-thresholded network."

th=0.826002

kinc run extract \
  --emx "PRJNA301554.slim.GEM.log2.emx" \
  --ccm "PRJNA301554.slim.GEM.log2.traditional.ccm" \
  --cmx "PRJNA301554.slim.GEM.log2.traditional.cmx" \
  --format "tidy" \
  --output "PRJNA301554.slim.GEM.log2.traditional.paf-th${th}-gcn.txt" \
  --mincorr  $th \
  --maxcorr 1
