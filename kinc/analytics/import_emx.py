"""
Import gene expression matrix from text file to HDF5.
"""

import numpy as np
import pandas as pd
from typing import Optional
from ..data import ExpressionMatrix


def import_emx(input_file: str, output_file: str,
               nan_token: str = 'NA',
               has_header: bool = True,
               has_row_names: bool = True) -> None:
    """
    Import gene expression matrix from text file.

    Args:
        input_file: Path to input text file (tab or space delimited)
        output_file: Path to output HDF5 file
        nan_token: String representing missing values
        has_header: Whether first row contains sample names
        has_row_names: Whether first column contains gene names
    """
    print(f"Importing expression matrix from {input_file}...")

    # Read text file with pandas
    if has_header and has_row_names:
        df = pd.read_csv(input_file, sep=r'\s+', index_col=0, na_values=[nan_token])
    elif has_header:
        df = pd.read_csv(input_file, sep=r'\s+', na_values=[nan_token])
    elif has_row_names:
        df = pd.read_csv(input_file, sep=r'\s+', index_col=0, header=None, na_values=[nan_token])
    else:
        df = pd.read_csv(input_file, sep=r'\s+', header=None, na_values=[nan_token])

    # Get dimensions
    n_genes, n_samples = df.shape

    # Get gene and sample names
    if has_row_names:
        gene_names = [str(name) for name in df.index]
    else:
        gene_names = [f"Gene_{i}" for i in range(n_genes)]

    if has_header:
        sample_names = [str(name) for name in df.columns]
    else:
        sample_names = [f"Sample_{i}" for i in range(n_samples)]

    print(f"  Genes: {n_genes}")
    print(f"  Samples: {n_samples}")

    # Convert to numpy array (float32)
    data = df.values.astype(np.float32)

    # Create expression matrix and write data
    emx = ExpressionMatrix(output_file, mode='w')
    emx.create(gene_names, sample_names)

    with emx:
        emx.write_expressions(data)

    print(f"Expression matrix saved to {output_file}")
