"""
CPMatrix - Cluster parameter matrix (GMM components) using HDF5.
"""

import numpy as np
import h5py
from typing import List, Optional


class CPMatrix:
    """
    Stores GMM cluster parameters in HDF5 format.

    Data structure:
        /pairs: group containing datasets indexed by pair_id
            pair_XXXXXXX/indices: (2,) int32 array [gene_x, gene_y]
            pair_XXXXXXX/components: structured array with fields:
                - pi: mixture weight (float32)
                - mu: mean vector (2, float32)
                - sigma: covariance matrix (2x2, float32)
        /genes: (n_genes,) string array
        /metadata: group with attributes
    """

    def __init__(self, filename: str, mode: str = 'r'):
        """Initialize CPMatrix."""
        self.filename = filename
        self.mode = mode
        self._file = None

    def __enter__(self):
        """Context manager entry."""
        self._file = h5py.File(self.filename, self.mode)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        if self._file is not None:
            self._file.close()
            self._file = None

    def create(self, gene_names: List[str]):
        """Create new CPMatrix file."""
        with h5py.File(self.filename, 'w') as f:
            f.create_group('pairs')
            f.create_dataset('genes',
                           data=np.array(gene_names, dtype=h5py.string_dtype()),
                           compression='gzip')
            metadata = f.create_group('metadata')
            metadata.attrs['version'] = '3.4.2'
            metadata.attrs['n_genes'] = len(gene_names)

    def __repr__(self):
        return f"CPMatrix(filename='{self.filename}')"
