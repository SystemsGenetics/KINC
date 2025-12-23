#!/bin/bash
# Run this after any KINC run to extract the various formats.

# Step 5: Extract the network.
echo "Extracting the condition-specific network."

p="1e-3"
r2="0.30"
th="0.00"

formats="tidy minimal graphml text"

for format in $formats; do

  kinc run extract \
    --emx "PRJNA301554.slim.GEM.log2.emx" \
    --ccm "PRJNA301554.slim.GEM.log2.paf.ccm" \
    --cmx "PRJNA301554.slim.GEM.log2.paf.cmx" \
    --csm "PRJNA301554.slim.GEM.log2.paf.csm" \
    --format $format \
    --output "PRJNA301554.slim.GEM.log2.paf-th${th}-p${p}-rsqr${r2}.${format}.txt" \
    --mincorr $th \
    --maxcorr 1 \
    --filter-pvalue $p \
    --filter-rsquare $r2

done;
