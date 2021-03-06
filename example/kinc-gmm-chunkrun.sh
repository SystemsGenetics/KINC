#!/bin/bash
# Perform a GMM run of KINC
mkdir -p results-kinc-gmm-chunkrun
cd results-kinc-gmm-chunkrun

# Initialize number of chunks
NP=4

# Initialize KINC settings
kinc settings set cuda 0
kinc settings set opencl none
kinc settings set threads 2
kinc settings set logging off

# Step 1: Import the expression matrix
echo "Importing the gene expression matrix (GEM) for the slimmed experiment PRJNA301554"

kinc run import-emx \
   --input "../data/PRJNA301554.slim.GEM.log2.txt" \
   --output "PRJNA301554.slim.GEM.log2.emx" \
   --nan "NA" \
   --samples 0

# Step 2: Create the similarity matrix.
echo "Performing similarity calculations using GMMs (GPUs make this go fastest)"

for i in $(seq 0 $(expr ${NP} - 1)); do
    echo "Running chunk ${i}"

    kinc chunkrun ${i} ${NP} similarity \
       --input "PRJNA301554.slim.GEM.log2.emx" \
       --ccm "PRJNA301554.slim.GEM.log2.ccm" \
       --cmx "PRJNA301554.slim.GEM.log2.cmx" \
       --clusmethod "gmm" \
       --corrmethod "spearman" \
       --minexpr -inf \
       --minsamp 25 \
       --minclus 1 \
       --maxclus 5 \
       --crit "ICL" \
       --preout TRUE \
       --postout TRUE \
       --mincorr 0 \
       --maxcorr 1
done

echo "Merging chunks"
kinc merge ${NP} similarity \
    --input "PRJNA301554.slim.GEM.log2.emx" \
    --ccm "PRJNA301554.slim.GEM.log2.ccm" \
    --cmx "PRJNA301554.slim.GEM.log2.cmx" \
    --clusmethod "gmm" \
    --corrmethod "spearman" \
    --minexpr -inf \
    --minsamp 25 \
    --minclus 1 \
    --maxclus 5 \
    --crit "ICL" \
    --preout TRUE \
    --postout TRUE \
    --mincorr 0 \
    --maxcorr 1

# Step 3: Filter for clusters with low power.
echo "Filtering edges with insufficient statistical power."

kinc run corrpower \
   --ccm-in "PRJNA301554.slim.GEM.log2.ccm" \
   --cmx-in "PRJNA301554.slim.GEM.log2.cmx" \
   --ccm-out "PRJNA301554.slim.GEM.log2.paf.ccm" \
   --cmx-out "PRJNA301554.slim.GEM.log2.paf.cmx" \
   --alpha 0.001 \
   --power 0.8

# Step 4: Condition-specific analysis.
echo "Performing condition-specific testing."

kinc run cond-test \
   --emx "PRJNA301554.slim.GEM.log2.emx" \
   --ccm "PRJNA301554.slim.GEM.log2.paf.ccm" \
   --cmx "PRJNA301554.slim.GEM.log2.paf.cmx" \
   --amx "../data/PRJNA301554.slim.annotations.txt" \
   --output "PRJNA301554.slim.GEM.log2.paf.csm" \
   --feat-tests "Subspecies,Treatment,GTAbbr,Duration" \
   --feat-types "Duration:quantitative" \
   --alpha 0.001 \
   --power 0.8

# Step 5: Extract the network.
echo "Extracting the condition-specific network."

p="1e-3"
r2="0.30"
th="0.00"

kinc run extract \
    --emx "PRJNA301554.slim.GEM.log2.emx" \
    --ccm "PRJNA301554.slim.GEM.log2.paf.ccm" \
    --cmx "PRJNA301554.slim.GEM.log2.paf.cmx" \
    --csm "PRJNA301554.slim.GEM.log2.paf.csm" \
    --format "tidy" \
    --output "PRJNA301554.slim.GEM.log2.paf-th${th}-p${p}-rsqr${r2}.txt" \
    --mincorr $th \
    --maxcorr 1 \
    --filter-pvalue $p \
    --filter-rsquare $r2

# Step 6: Removed Biased Edges
echo "Removing biased edges from the extract network using KINC.R."

kinc-filter-bias.R \
    --net "PRJNA301554.slim.GEM.log2.paf-th0.00-p1e-3-rsqr0.30.txt" \
    --emx "../data/PRJNA301554.slim.GEM.log2.txt" \
    --out_prefix "PRJNA301554.slim.GEM.log2.paf-th0.00-p1e-3-rsqr0.30"

# Step 7: Generate summary plots
echo "Generating summary plots using KINC.R."

kinc-make-plots.R \
    --net "PRJNA301554.slim.GEM.log2.paf-th0.00-p1e-3-rsqr0.30-filtered.GCN.txt" \
    --out_prefix "PRJNA301554.slim.GEM.log2.paf-th0.00-p1e-3-rsqr0.30-filtered"

# Step 8: Threshold the network by ranks and generate condition-specific networks
echo "Generating condition-specific network files by class and label."

kinc-filter-rank.R \
    --net "PRJNA301554.slim.GEM.log2.paf-th0.00-p1e-3-rsqr0.30-filtered.GCN.txt" \
    --out_prefix "PRJNA301554.slim.GEM.log2.paf-th0.00-p1e-3-rsqr0.30-filtered" \
    --top_n 25000

kinc-filter-rank.R \
    --net "PRJNA301554.slim.GEM.log2.paf-th0.00-p1e-3-rsqr0.30-filtered.GCN.txt" \
    --out_prefix "PRJNA301554.slim.GEM.log2.paf-th0.00-p1e-3-rsqr0.30-filtered" \
    --save_condition_networks \
    --top_n 25000

kinc-filter-rank.R \
    --net "PRJNA301554.slim.GEM.log2.paf-th0.00-p1e-3-rsqr0.30-filtered.GCN.txt" \
    --out_prefix "PRJNA301554.slim.GEM.log2.paf-th0.00-p1e-3-rsqr0.30-filtered" \
    --save_condition_networks --unique_filter "label" \
    --top_n 25000

kinc-filter-rank.R \
    --net "PRJNA301554.slim.GEM.log2.paf-th0.00-p1e-3-rsqr0.30-filtered.GCN.txt" \
    --out_prefix "PRJNA301554.slim.GEM.log2.paf-th0.00-p1e-3-rsqr0.30-filtered" \
    --save_condition_networks --unique_filter "class" \
    --top_n 25000
