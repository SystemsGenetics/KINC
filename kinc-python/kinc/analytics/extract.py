"""
Extract network edges from correlation matrix to text file.
"""

import numpy as np
from typing import Optional
from ..data import CorrelationMatrix, CCMatrix


def extract_network(cmx_file: str,
                    output_file: str,
                    ccm_file: Optional[str] = None,
                    min_correlation: float = 0.0,
                    max_correlation: float = 1.0,
                    output_format: str = 'text') -> None:
    """
    Extract network edges from correlation matrix.

    Args:
        cmx_file: Input correlation matrix (HDF5)
        output_file: Output network file
        ccm_file: Optional cluster composition matrix (HDF5)
        min_correlation: Minimum |correlation| threshold
        max_correlation: Maximum |correlation| threshold
        output_format: 'text', 'tidy', or 'minimal'
    """
    print(f"Extracting network from {cmx_file}...")

    # Load correlation matrix
    cmx = CorrelationMatrix(cmx_file, mode='r')
    gene_names = cmx.gene_names

    # Optionally load CCMatrix
    if ccm_file:
        ccm = CCMatrix(ccm_file, mode='r')
        sample_names = ccm.sample_names
    else:
        ccm = None
        sample_names = None

    # Open output file
    with open(output_file, 'w') as f:
        # Write header based on format
        if output_format == 'text':
            if ccm:
                f.write("Source\tTarget\tCluster\tCorrelation\tSampleMask\n")
            else:
                f.write("Source\tTarget\tCluster\tCorrelation\n")
        elif output_format == 'tidy':
            f.write("Gene1\tGene2\tCluster\tCorrelation\n")
        else:  # minimal
            f.write("Gene1\tGene2\tCorrelation\n")

        # Iterate over all pairs
        edge_count = 0
        for gene_x, gene_y, correlations in cmx.iter_pairs():
            gene_x_name = gene_names[gene_x]
            gene_y_name = gene_names[gene_y]

            # Get sample masks if available
            if ccm:
                sample_masks = ccm.read_pair(gene_x, gene_y)
            else:
                sample_masks = None

            # Write edges for each cluster
            for cluster_idx, corr in enumerate(correlations):
                # Check threshold
                if np.isnan(corr):
                    continue

                abs_corr = abs(corr)
                if abs_corr < min_correlation or abs_corr > max_correlation:
                    continue

                # Format output based on format type
                if output_format == 'text':
                    if sample_masks is not None:
                        mask_str = ''.join(['1' if m == cluster_idx else '0'
                                          for m in sample_masks[cluster_idx]])
                        f.write(f"{gene_x_name}\t{gene_y_name}\t{cluster_idx}\t{corr:.6f}\t{mask_str}\n")
                    else:
                        f.write(f"{gene_x_name}\t{gene_y_name}\t{cluster_idx}\t{corr:.6f}\n")
                elif output_format == 'tidy':
                    f.write(f"{gene_x_name}\t{gene_y_name}\t{cluster_idx}\t{corr:.6f}\n")
                else:  # minimal
                    if cluster_idx == 0:  # Only write first cluster for minimal format
                        f.write(f"{gene_x_name}\t{gene_y_name}\t{corr:.6f}\n")

                edge_count += 1

    print(f"Network extracted to {output_file}")
    print(f"  Total edges: {edge_count:,}")
