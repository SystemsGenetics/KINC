"""
ExpressionMatrix - Gene expression data storage using HDF5.
"""

import numpy as np
import h5py
from typing import Optional, List


class ExpressionMatrix:
    """
    Stores gene expression data in HDF5 format.

    Data structure:
        /expressions: (n_genes, n_samples) float32 array
        /genes: (n_genes,) string array of gene names
        /samples: (n_samples,) string array of sample names
        /metadata: group containing additional attributes

    Attributes:
        filename: Path to HDF5 file
        n_genes: Number of genes
        n_samples: Number of samples
        gene_names: List of gene names
        sample_names: List of sample names
    """

    def __init__(self, filename: str, mode: str = 'r'):
        """
        Initialize ExpressionMatrix.

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
            self._n_genes = f['expressions'].shape[0]
            self._n_samples = f['expressions'].shape[1]
            self._gene_names = [g.decode('utf-8') if isinstance(g, bytes) else g
                               for g in f['genes'][:]]
            self._sample_names = [s.decode('utf-8') if isinstance(s, bytes) else s
                                 for s in f['samples'][:]]

    def create(self, gene_names: List[str], sample_names: List[str]):
        """
        Create new expression matrix file.

        Args:
            gene_names: List of gene names
            sample_names: List of sample names
        """
        self._n_genes = len(gene_names)
        self._n_samples = len(sample_names)
        self._gene_names = gene_names
        self._sample_names = sample_names

        with h5py.File(self.filename, 'w') as f:
            # Create datasets
            f.create_dataset('expressions',
                           shape=(self._n_genes, self._n_samples),
                           dtype='float32',
                           chunks=(min(1000, self._n_genes), self._n_samples),
                           compression='gzip',
                           compression_opts=4)

            # Store gene and sample names
            f.create_dataset('genes',
                           data=np.array(gene_names, dtype=h5py.string_dtype()),
                           compression='gzip')
            f.create_dataset('samples',
                           data=np.array(sample_names, dtype=h5py.string_dtype()),
                           compression='gzip')

            # Create metadata group
            metadata = f.create_group('metadata')
            metadata.attrs['version'] = '3.4.2'
            metadata.attrs['n_genes'] = self._n_genes
            metadata.attrs['n_samples'] = self._n_samples

    def write_expression(self, gene_idx: int, values: np.ndarray):
        """
        Write expression values for a single gene.

        Args:
            gene_idx: Gene index
            values: Expression values (n_samples,)
        """
        if self._file is None:
            raise RuntimeError("File not opened. Use 'with' statement or open() method.")

        self._file['expressions'][gene_idx, :] = values

    def write_expressions(self, data: np.ndarray):
        """
        Write entire expression matrix at once.

        Args:
            data: Expression matrix (n_genes, n_samples)
        """
        if self._file is None:
            raise RuntimeError("File not opened. Use 'with' statement or open() method.")

        self._file['expressions'][:, :] = data

    def read_expression(self, gene_idx: int) -> np.ndarray:
        """
        Read expression values for a single gene.

        Args:
            gene_idx: Gene index

        Returns:
            Expression values (n_samples,)
        """
        if self._file is None:
            with h5py.File(self.filename, 'r') as f:
                return f['expressions'][gene_idx, :].astype(np.float32)
        else:
            return self._file['expressions'][gene_idx, :].astype(np.float32)

    def read_expressions(self, gene_indices: Optional[List[int]] = None) -> np.ndarray:
        """
        Read expression values for multiple genes.

        Args:
            gene_indices: List of gene indices (None = all genes)

        Returns:
            Expression matrix (n_genes, n_samples) or subset
        """
        if self._file is None:
            with h5py.File(self.filename, 'r') as f:
                if gene_indices is None:
                    return f['expressions'][:, :].astype(np.float32)
                else:
                    return f['expressions'][gene_indices, :].astype(np.float32)
        else:
            if gene_indices is None:
                return self._file['expressions'][:, :].astype(np.float32)
            else:
                return self._file['expressions'][gene_indices, :].astype(np.float32)

    def read_pair(self, gene_x: int, gene_y: int) -> tuple:
        """
        Read expression values for a gene pair.

        Args:
            gene_x: First gene index
            gene_y: Second gene index

        Returns:
            Tuple of (expr_x, expr_y) arrays
        """
        if self._file is None:
            with h5py.File(self.filename, 'r') as f:
                expr_x = f['expressions'][gene_x, :].astype(np.float32)
                expr_y = f['expressions'][gene_y, :].astype(np.float32)
                return expr_x, expr_y
        else:
            expr_x = self._file['expressions'][gene_x, :].astype(np.float32)
            expr_y = self._file['expressions'][gene_y, :].astype(np.float32)
            return expr_x, expr_y

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

    @property
    def gene_names(self) -> List[str]:
        """Gene names."""
        if self._gene_names is None:
            self._load_metadata()
        return self._gene_names

    @property
    def sample_names(self) -> List[str]:
        """Sample names."""
        if self._sample_names is None:
            self._load_metadata()
        return self._sample_names

    def __repr__(self):
        return f"ExpressionMatrix(filename='{self.filename}', n_genes={self.n_genes}, n_samples={self.n_samples})"
