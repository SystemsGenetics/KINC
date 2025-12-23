"""
CorrelationMatrix - Pairwise correlation storage using HDF5.
"""

import numpy as np
import h5py
from typing import Optional, List, Tuple


class CorrelationMatrix:
    """
    Stores pairwise gene correlations in HDF5 format with support for multiple
    clusters per gene pair (from GMM clustering).

    Data structure:
        /pairs: group containing datasets indexed by pair_id
            pair_XXXXXXX/indices: (2,) int32 array [gene_x, gene_y]
            pair_XXXXXXX/correlations: (n_clusters,) float32 array
        /genes: (n_genes,) string array of gene names
        /metadata: group containing attributes
            - correlation_method: 'pearson' or 'spearman'
            - n_genes: number of genes
            - max_clusters: maximum clusters per pair

    Attributes:
        filename: Path to HDF5 file
        n_genes: Number of genes
        correlation_method: Correlation method used
    """

    def __init__(self, filename: str, mode: str = 'r'):
        """
        Initialize CorrelationMatrix.

        Args:
            filename: Path to HDF5 file
            mode: File mode ('r', 'w', 'a')
        """
        self.filename = filename
        self.mode = mode
        self._file = None
        self._n_genes = None
        self._gene_names = None
        self._correlation_method = None
        self._max_clusters = None
        self._pair_count = 0

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
            self._correlation_method = metadata.attrs['correlation_method']
            self._max_clusters = metadata.attrs.get('max_clusters', 5)
            self._gene_names = [g.decode('utf-8') if isinstance(g, bytes) else g
                               for g in f['genes'][:]]
            self._pair_count = len(f['pairs'].keys()) if 'pairs' in f else 0

    def create(self, gene_names: List[str], correlation_method: str = 'pearson',
               max_clusters: int = 5):
        """
        Create new correlation matrix file.

        Args:
            gene_names: List of gene names
            correlation_method: 'pearson' or 'spearman'
            max_clusters: Maximum clusters per pair
        """
        self._n_genes = len(gene_names)
        self._gene_names = gene_names
        self._correlation_method = correlation_method
        self._max_clusters = max_clusters

        with h5py.File(self.filename, 'w') as f:
            # Create pairs group
            f.create_group('pairs')

            # Store gene names
            f.create_dataset('genes',
                           data=np.array(gene_names, dtype=h5py.string_dtype()),
                           compression='gzip')

            # Create metadata
            metadata = f.create_group('metadata')
            metadata.attrs['version'] = '3.4.2'
            metadata.attrs['n_genes'] = self._n_genes
            metadata.attrs['correlation_method'] = correlation_method
            metadata.attrs['max_clusters'] = max_clusters

    def write_pair(self, gene_x: int, gene_y: int, correlations: np.ndarray):
        """
        Write correlation values for a gene pair.

        Args:
            gene_x: First gene index (must be > gene_y)
            gene_y: Second gene index
            correlations: Correlation values for each cluster
        """
        if self._file is None:
            raise RuntimeError("File not opened. Use 'with' statement or open() method.")

        if gene_x <= gene_y:
            raise ValueError("gene_x must be > gene_y (lower triangle only)")

        # Create pair identifier (lower triangle index)
        pair_id = self._get_pair_id(gene_x, gene_y)
        pair_key = f"pair_{pair_id:08d}"

        # Create group for this pair if it doesn't exist
        if pair_key not in self._file['pairs']:
            pair_group = self._file['pairs'].create_group(pair_key)

            # Store indices
            pair_group.create_dataset('indices',
                                    data=np.array([gene_x, gene_y], dtype=np.int32))

            # Store correlations
            pair_group.create_dataset('correlations',
                                    data=correlations.astype(np.float32),
                                    compression='gzip')
            self._pair_count += 1

    def read_pair(self, gene_x: int, gene_y: int) -> Optional[np.ndarray]:
        """
        Read correlation values for a gene pair.

        Args:
            gene_x: First gene index
            gene_y: Second gene index

        Returns:
            Correlation values array or None if pair not found
        """
        # Ensure lower triangle
        if gene_x < gene_y:
            gene_x, gene_y = gene_y, gene_x

        pair_id = self._get_pair_id(gene_x, gene_y)
        pair_key = f"pair_{pair_id:08d}"

        if self._file is None:
            with h5py.File(self.filename, 'r') as f:
                if pair_key in f['pairs']:
                    return f['pairs'][pair_key]['correlations'][:].astype(np.float32)
                return None
        else:
            if pair_key in self._file['pairs']:
                return self._file['pairs'][pair_key]['correlations'][:].astype(np.float32)
            return None

    def iter_pairs(self):
        """
        Iterate over all gene pairs.

        Yields:
            Tuple of (gene_x, gene_y, correlations)
        """
        with h5py.File(self.filename, 'r') as f:
            for pair_key in f['pairs'].keys():
                pair_group = f['pairs'][pair_key]
                indices = pair_group['indices'][:]
                correlations = pair_group['correlations'][:].astype(np.float32)
                yield int(indices[0]), int(indices[1]), correlations

    def _get_pair_id(self, gene_x: int, gene_y: int) -> int:
        """
        Calculate unique pair ID from gene indices (lower triangle only).

        Args:
            gene_x: First gene index (x > y)
            gene_y: Second gene index

        Returns:
            Unique pair ID
        """
        return gene_x * (gene_x - 1) // 2 + gene_y

    @property
    def n_genes(self) -> int:
        """Number of genes."""
        if self._n_genes is None:
            self._load_metadata()
        return self._n_genes

    @property
    def correlation_method(self) -> str:
        """Correlation method."""
        if self._correlation_method is None:
            self._load_metadata()
        return self._correlation_method

    @property
    def gene_names(self) -> List[str]:
        """Gene names."""
        if self._gene_names is None:
            self._load_metadata()
        return self._gene_names

    @property
    def pair_count(self) -> int:
        """Number of pairs stored."""
        if self._file is None:
            with h5py.File(self.filename, 'r') as f:
                return len(f['pairs'].keys())
        else:
            return len(self._file['pairs'].keys())

    def __repr__(self):
        return f"CorrelationMatrix(filename='{self.filename}', n_genes={self.n_genes}, pairs={self.pair_count}, method='{self.correlation_method}')"
