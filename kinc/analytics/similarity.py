"""
Similarity analytic - Pairwise gene correlation with optional GMM clustering.
"""

import numpy as np
from tqdm import tqdm
from typing import Optional, Tuple, List
from ..data import ExpressionMatrix, CorrelationMatrix, CCMatrix
from ..core import GMMClusterer, pearson_correlation, spearman_correlation, compute_correlations
from ..core.outliers import remove_outliers_all_clusters


class Similarity:
    """
    Compute pairwise gene correlations with optional GMM clustering.

    This is the core KINC analytic that computes similarity (correlation)
    between all gene pairs, with optional Gaussian Mixture Model clustering
    to identify condition-specific correlations.
    """

    def __init__(self,
                 input_file: str,
                 cmx_file: str,
                 ccm_file: str,
                 correlation_method: str = 'pearson',
                 clustering_method: str = 'none',
                 min_expression: float = -np.inf,
                 max_expression: float = np.inf,
                 min_samples: int = 30,
                 min_clusters: int = 1,
                 max_clusters: int = 5,
                 criterion: str = 'ICL',
                 preout: bool = True,
                 postout: bool = True,
                 min_correlation: float = 0.5,
                 max_correlation: float = 1.0):
        """
        Initialize Similarity analytic.

        Args:
            input_file: Input expression matrix (HDF5)
            cmx_file: Output correlation matrix (HDF5)
            ccm_file: Output cluster composition matrix (HDF5)
            correlation_method: 'pearson' or 'spearman'
            clustering_method: 'none' or 'gmm'
            min_expression: Minimum expression threshold
            max_expression: Maximum expression threshold
            min_samples: Minimum samples required for correlation
            min_clusters: Minimum clusters to test (GMM)
            max_clusters: Maximum clusters to test (GMM)
            criterion: Model selection criterion ('AIC', 'BIC', 'ICL')
            preout: Remove pre-clustering outliers
            postout: Remove post-clustering outliers
            min_correlation: Minimum |correlation| to save
            max_correlation: Maximum |correlation| to save
        """
        self.input_file = input_file
        self.cmx_file = cmx_file
        self.ccm_file = ccm_file
        self.correlation_method = correlation_method.lower()
        self.clustering_method = clustering_method.lower()
        self.min_expression = min_expression
        self.max_expression = max_expression
        self.min_samples = min_samples
        self.min_clusters = min_clusters
        self.max_clusters = max_clusters
        self.criterion = criterion
        self.preout = preout
        self.postout = postout
        self.min_correlation = min_correlation
        self.max_correlation = max_correlation

        # Validate parameters
        if self.correlation_method not in ['pearson', 'spearman']:
            raise ValueError(f"Invalid correlation method: {correlation_method}")
        if self.clustering_method not in ['none', 'gmm']:
            raise ValueError(f"Invalid clustering method: {clustering_method}")

        # Load expression matrix metadata
        self.emx = ExpressionMatrix(input_file, mode='r')
        self.n_genes = self.emx.n_genes
        self.n_samples = self.emx.n_samples
        self.gene_names = self.emx.gene_names
        self.sample_names = self.emx.sample_names

        # Calculate total pairs
        self.total_pairs = self.n_genes * (self.n_genes - 1) // 2

        print(f"Similarity Configuration:")
        print(f"  Genes: {self.n_genes}")
        print(f"  Samples: {self.n_samples}")
        print(f"  Total pairs: {self.total_pairs:,}")
        print(f"  Correlation method: {self.correlation_method}")
        print(f"  Clustering method: {self.clustering_method}")
        if self.clustering_method == 'gmm':
            print(f"  Cluster range: {self.min_clusters}-{self.max_clusters}")
            print(f"  Criterion: {self.criterion}")
        print(f"  Min samples: {self.min_samples}")
        print(f"  Outlier removal: pre={self.preout}, post={self.postout}")
        print(f"  Correlation filter: [{self.min_correlation}, {self.max_correlation}]")

    def run(self, n_jobs: int = 1) -> None:
        """
        Run similarity computation.

        Args:
            n_jobs: Number of parallel jobs (1 = serial, >1 = parallel)
        """
        # Initialize output matrices
        self._initialize_outputs()

        if n_jobs == 1:
            # Serial execution
            self._run_serial()
        else:
            # Parallel execution
            self._run_parallel(n_jobs)

        print(f"\nSimilarity computation complete!")
        print(f"  Pairs processed: {self.total_pairs:,}")
        print(f"  Pairs saved: {self.pairs_saved:,}")

    def _initialize_outputs(self):
        """Initialize output correlation and CCMatrix files."""
        # Create correlation matrix
        cmx = CorrelationMatrix(self.cmx_file, mode='w')
        cmx.create(self.gene_names, self.correlation_method, self.max_clusters)

        # Create CCMatrix
        ccm = CCMatrix(self.ccm_file, mode='w')
        ccm.create(self.gene_names, self.sample_names)

        self.pairs_saved = 0

    def _run_serial(self):
        """Run similarity computation in serial mode."""
        print("\nRunning in serial mode...")

        # Open output files
        with CorrelationMatrix(self.cmx_file, mode='a') as cmx, \
             CCMatrix(self.ccm_file, mode='a') as ccm, \
             ExpressionMatrix(self.input_file, mode='r') as emx:

            # Iterate over all gene pairs
            pair_count = 0
            with tqdm(total=self.total_pairs, desc="Processing pairs") as pbar:
                for gene_x in range(1, self.n_genes):
                    for gene_y in range(gene_x):
                        # Process this pair
                        result = self._process_pair(emx, gene_x, gene_y)

                        if result is not None:
                            n_clusters, correlations, sample_masks = result

                            # Save if correlations meet threshold
                            if self._should_save(correlations):
                                cmx.write_pair(gene_x, gene_y, correlations)
                                ccm.write_pair(gene_x, gene_y, sample_masks)
                                self.pairs_saved += 1

                        pair_count += 1
                        pbar.update(1)

    def _process_pair(self, emx: ExpressionMatrix,
                     gene_x: int, gene_y: int) -> Optional[Tuple[int, np.ndarray, np.ndarray]]:
        """
        Process a single gene pair.

        Args:
            emx: Expression matrix
            gene_x: First gene index
            gene_y: Second gene index

        Returns:
            Tuple of (n_clusters, correlations, sample_masks) or None if insufficient data
        """
        # Read expression data for this pair
        expr_x, expr_y = emx.read_pair(gene_x, gene_y)

        # Initialize sample labels
        # 0 = valid, -6 = threshold, -7 = pre-outlier, -8 = post-outlier, -9 = NaN
        labels = np.zeros(self.n_samples, dtype=np.int8)

        # Mark NaN values
        nan_mask = np.isnan(expr_x) | np.isnan(expr_y)
        labels[nan_mask] = -9

        # Mark expression threshold violations
        threshold_mask = ((expr_x < self.min_expression) | (expr_x > self.max_expression) |
                         (expr_y < self.min_expression) | (expr_y > self.max_expression))
        labels[threshold_mask & ~nan_mask] = -6

        # Valid samples for clustering
        valid_mask = labels == 0

        # Check minimum samples
        if np.sum(valid_mask) < self.min_samples:
            return None

        # Pre-clustering outlier removal
        if self.preout:
            labels = self._remove_outliers_single_cluster(expr_x, expr_y, labels)
            valid_mask = labels == 0
            if np.sum(valid_mask) < self.min_samples:
                return None

        # Clustering
        if self.clustering_method == 'gmm':
            n_clusters, labels = self._cluster_gmm(expr_x, expr_y, valid_mask, labels)
        else:
            n_clusters = 1
            # Labels already initialized to 0 for valid samples

        # Post-clustering outlier removal
        if self.postout and self.clustering_method == 'gmm':
            labels = remove_outliers_all_clusters(expr_x, expr_y, labels, n_clusters)

        # Compute correlations for each cluster
        correlations = compute_correlations(
            expr_x, expr_y, labels, n_clusters,
            method=self.correlation_method,
            min_samples=self.min_samples
        )

        # Create sample masks (one row per cluster)
        sample_masks = np.zeros((n_clusters, self.n_samples), dtype=np.int8)
        for k in range(n_clusters):
            sample_masks[k, :] = (labels == k).astype(np.int8)

        return n_clusters, correlations, sample_masks

    def _cluster_gmm(self, expr_x: np.ndarray, expr_y: np.ndarray,
                     valid_mask: np.ndarray, labels: np.ndarray) -> Tuple[int, np.ndarray]:
        """Apply GMM clustering."""
        clusterer = GMMClusterer(
            min_clusters=self.min_clusters,
            max_clusters=self.max_clusters,
            criterion=self.criterion
        )

        n_clusters, new_labels = clusterer.fit(expr_x, expr_y, valid_mask)

        # Merge new cluster labels with existing labels
        # (preserving -6, -7, -8, -9 markers)
        final_labels = labels.copy()
        final_labels[valid_mask] = new_labels[valid_mask]

        return n_clusters, final_labels

    def _remove_outliers_single_cluster(self, expr_x: np.ndarray,
                                       expr_y: np.ndarray,
                                       labels: np.ndarray) -> np.ndarray:
        """Remove outliers from cluster 0 (all valid samples)."""
        from ..core.outliers import remove_outliers
        return remove_outliers(expr_x, expr_y, labels, cluster=0)

    def _should_save(self, correlations: np.ndarray) -> bool:
        """Check if any correlation meets the threshold criteria."""
        abs_corrs = np.abs(correlations)
        valid_corrs = abs_corrs[~np.isnan(abs_corrs)]

        if len(valid_corrs) == 0:
            return False

        # Check if any correlation is in the acceptable range
        in_range = (valid_corrs >= self.min_correlation) & (valid_corrs <= self.max_correlation)
        return np.any(in_range)

    def _run_parallel(self, n_jobs: int):
        """Run similarity computation in parallel mode (TODO)."""
        raise NotImplementedError("Parallel execution not yet implemented. Use n_jobs=1 for serial mode.")
