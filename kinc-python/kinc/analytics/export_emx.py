"""
Export gene expression matrix from HDF5 to text file.
"""

import numpy as np
import pandas as pd
from ..data import ExpressionMatrix


def export_emx(input_file: str, output_file: str,
               nan_token: str = 'NA') -> None:
    """
    Export gene expression matrix to text file.

    Args:
        input_file: Path to input HDF5 file
        output_file: Path to output text file
        nan_token: String to use for missing values
    """
    print(f"Exporting expression matrix from {input_file}...")

    # Load expression matrix
    emx = ExpressionMatrix(input_file, mode='r')

    # Read all data
    data = emx.read_expressions()
    gene_names = emx.gene_names
    sample_names = emx.sample_names

    # Create DataFrame
    df = pd.DataFrame(data, index=gene_names, columns=sample_names)

    # Write to text file
    df.to_csv(output_file, sep='\t', na_rep=nan_token, float_format='%.6f')

    print(f"Expression matrix exported to {output_file}")
    print(f"  Genes: {len(gene_names)}")
    print(f"  Samples: {len(sample_names)}")
