#!/usr/bin/env python
"""
Example KINC-Python workflow using Python API.
"""

from kinc.analytics import import_emx, export_emx, Similarity, extract_network

# Input/output files
INPUT_GEM = "data.txt"
EMX_FILE = "data.emx.h5"
CMX_FILE = "data.cmx.h5"
CCM_FILE = "data.ccm.h5"
NETWORK_FILE = "network.txt"

print("=== KINC-Python Example Workflow ===\n")

# Step 1: Import expression matrix
print("Step 1: Importing expression matrix...")
import_emx(
    input_file=INPUT_GEM,
    output_file=EMX_FILE,
    nan_token="NA"
)

# Step 2: Compute similarity with GMM clustering
print("\nStep 2: Computing similarity with GMM clustering...")
sim = Similarity(
    input_file=EMX_FILE,
    cmx_file=CMX_FILE,
    ccm_file=CCM_FILE,
    correlation_method='spearman',
    clustering_method='gmm',
    min_clusters=1,
    max_clusters=5,
    criterion='ICL',
    min_samples=30,
    preout=True,
    postout=True,
    min_correlation=0.85
)
sim.run(n_jobs=1)

# Step 3: Extract network
print("\nStep 3: Extracting network edges...")
extract_network(
    cmx_file=CMX_FILE,
    output_file=NETWORK_FILE,
    ccm_file=CCM_FILE,
    min_correlation=0.85,
    output_format='text'
)

print("\n=== Workflow complete! ===")
print(f"Network saved to: {NETWORK_FILE}")
