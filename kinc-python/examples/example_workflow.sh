#!/bin/bash
#
# Example KINC-Python workflow
#
# This script demonstrates a complete KINC workflow:
# 1. Import expression matrix from text
# 2. Compute similarity with GMM clustering
# 3. Extract network edges
#

# Input/output files
INPUT_GEM="data.txt"
EMX_FILE="data.emx.h5"
CMX_FILE="data.cmx.h5"
CCM_FILE="data.ccm.h5"
NETWORK_FILE="network.txt"

echo "=== KINC-Python Example Workflow ==="
echo

# Step 1: Import expression matrix
echo "Step 1: Importing expression matrix..."
kinc import-emx \
  --input $INPUT_GEM \
  --output $EMX_FILE \
  --nan "NA"

echo
echo "Step 2: Computing similarity with GMM clustering..."
kinc similarity \
  --input $EMX_FILE \
  --cmx $CMX_FILE \
  --ccm $CCM_FILE \
  --corrmethod spearman \
  --clusmethod gmm \
  --minclus 1 \
  --maxclus 5 \
  --crit ICL \
  --minsamp 30 \
  --preout true \
  --postout true \
  --mincorr 0.85

echo
echo "Step 3: Extracting network edges..."
kinc extract \
  --input $CMX_FILE \
  --ccm $CCM_FILE \
  --output $NETWORK_FILE \
  --mincorr 0.85 \
  --format text

echo
echo "=== Workflow complete! ==="
echo "Network saved to: $NETWORK_FILE"
