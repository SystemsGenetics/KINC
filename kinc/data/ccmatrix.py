"""
CCMatrix - Cluster composition matrix (sample masks) using HDF5.
"""

import numpy as np
import h5py
from typing import List, Optional


class CCMatrix:
    """
    Stores cluster sample masks in HDF5 format.

    Data structure:
        /pairs: group containing datasets indexed by pair_id
            pair_XXXXXXX/indices: (2,) int32 array [gene_x, gene_y]
            pair_XXXXXXX/sample_masks: (n_clusters, n_samples) int8 array
        /genes: (n_genes,) string array
        /samples: (n_samples,) string array
        /metadata: group with attributes

    Sample mask values:
        0-15: Cluster assignment (0 = cluster 0, 1 = cluster 1, etc.)
        -6: Below/above expression threshold
        -7: Pre-clustering outlier
        -8: Post-clustering outlier
        -9: Missing value (NaN)
    """

    def __init__(self, filename: str, mode: str = 'r'):
        """
        Initialize CCMatrix.

        Args:
            filename: Path to HDF5 file
            mode: File mode ('r', 'w', 'a')
        """
        self.filename = filename
        self.mode = mode
        self._file = None
        self._n_genes = None
        self._n_samples = None
        self._gene_names = None
        self._sample_names = None

        if mode == 'r':
            self._load_metadata()

    def __enter__(self):
        """Context manager entry."""
        self._file = h5py.File(self.filename, self.mode)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        if self._file is not None:
            self._file.close()
            self._file = None

    def _load_metadata(self):
        """Load metadata from existing file."""
        with h5py.File(self.filename, 'r') as f:
            metadata = f['metadata']
            self._n_genes = metadata.attrs['n_genes']
            self._n_samples = metadata.attrs['n_samples']
            self._gene_names = [g.decode('utf-8') if isinstance(g, bytes) else g
                               for g in f['genes'][:]]
            self._sample_names = [s.decode('utf-8') if isinstance(s, bytes) else s
                                 for s in f['samples'][:]]

    def create(self, gene_names: List[str], sample_names: List[str]):
        """
        Create new CCMatrix file.

        Args:
            gene_names: List of gene names
            sample_names: List of sample names
        """
        self._n_genes = len(gene_names)
        self._n_samples = len(sample_names)
        self._gene_names = gene_names
        self._sample_names = sample_names

        with h5py.File(self.filename, 'w') as f:
            # Create pairs group
            f.create_group('pairs')

            # Store names
            f.create_dataset('genes',
                           data=np.array(gene_names, dtype=h5py.string_dtype()),
                           compression='gzip')
            f.create_dataset('samples',
                           data=np.array(sample_names, dtype=h5py.string_dtype()),
                           compression='gzip')

            # Create metadata
            metadata = f.create_group('metadata')
            metadata.attrs['version'] = '3.4.2'
            metadata.attrs['n_genes'] = self._n_genes
            metadata.attrs['n_samples'] = self._n_samples

    def write_pair(self, gene_x: int, gene_y: int, sample_masks: np.ndarray):
        """
        Write sample masks for a gene pair.

        Args:
            gene_x: First gene index (must be > gene_y)
            gene_y: Second gene index
            sample_masks: Sample masks array (n_clusters, n_samples) int8
        """
        if self._file is None:
            raise RuntimeError("File not opened. Use 'with' statement or open() method.")

        if gene_x <= gene_y:
            raise ValueError("gene_x must be > gene_y (lower triangle only)")

        pair_id = gene_x * (gene_x - 1) // 2 + gene_y
        pair_key = f"pair_{pair_id:08d}"

        if pair_key not in self._file['pairs']:
            pair_group = self._file['pairs'].create_group(pair_key)
            pair_group.create_dataset('indices',
                                    data=np.array([gene_x, gene_y], dtype=np.int32))
            pair_group.create_dataset('sample_masks',
                                    data=sample_masks.astype(np.int8),
                                    compression='gzip')

    def read_pair(self, gene_x: int, gene_y: int) -> Optional[np.ndarray]:
        """
        Read sample masks for a gene pair.

        Args:
            gene_x: First gene index
            gene_y: Second gene index

        Returns:
            Sample masks array (n_clusters, n_samples) or None
        """
        if gene_x < gene_y:
            gene_x, gene_y = gene_y, gene_x

        pair_id = gene_x * (gene_x - 1) // 2 + gene_y
        pair_key = f"pair_{pair_id:08d}"

        if self._file is None:
            with h5py.File(self.filename, 'r') as f:
                if pair_key in f['pairs']:
                    return f['pairs'][pair_key]['sample_masks'][:].astype(np.int8)
                return None
        else:
            if pair_key in self._file['pairs']:
                return self._file['pairs'][pair_key]['sample_masks'][:].astype(np.int8)
            return None

    @property
    def n_genes(self) -> int:
        """Number of genes."""
        if self._n_genes is None:
            self._load_metadata()
        return self._n_genes

    @property
    def n_samples(self) -> int:
        """Number of samples."""
        if self._n_samples is None:
            self._load_metadata()
        return self._n_samples

    def __repr__(self):
        return f"CCMatrix(filename='{self.filename}', n_genes={self.n_genes}, n_samples={self.n_samples})"
