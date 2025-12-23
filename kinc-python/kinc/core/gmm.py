"""
Gaussian Mixture Model clustering for gene pairs.
"""

import numpy as np
from sklearn.mixture import GaussianMixture
from typing import Tuple


class GMMClusterer:
    """
    GMM-based clustering for gene pairs with model selection.

    Uses scikit-learn's GaussianMixture with custom model selection
    supporting AIC, BIC, and ICL criteria.
    """

    def __init__(self, min_clusters: int = 1, max_clusters: int = 5,
                 criterion: str = 'ICL', random_state: int = 0):
        """
        Initialize GMM clusterer.

        Args:
            min_clusters: Minimum number of clusters to test
            max_clusters: Maximum number of clusters to test
            criterion: Model selection criterion ('AIC', 'BIC', 'ICL')
            random_state: Random seed for reproducibility
        """
        self.min_clusters = min_clusters
        self.max_clusters = max_clusters
        self.criterion = criterion.upper()
        self.random_state = random_state

        if self.criterion not in ['AIC', 'BIC', 'ICL']:
            raise ValueError(f"Invalid criterion: {criterion}. Must be AIC, BIC, or ICL.")

    def fit(self, x: np.ndarray, y: np.ndarray,
            sample_mask: np.ndarray) -> Tuple[int, np.ndarray]:
        """
        Fit GMM to gene pair data and select optimal number of clusters.

        Args:
            x: Expression values for gene X
            y: Expression values for gene Y
            sample_mask: Boolean mask of valid samples

        Returns:
            Tuple of (n_clusters, labels)
                n_clusters: Optimal number of clusters
                labels: Cluster assignments for each sample
        """
        # Extract valid samples
        x_valid = x[sample_mask]
        y_valid = y[sample_mask]

        # Stack into 2D array (n_samples, 2)
        X = np.column_stack([x_valid, y_valid])

        # Check if we have enough samples
        n_samples = X.shape[0]
        if n_samples < self.min_clusters * 3:  # Heuristic: at least 3 samples per cluster
            # Not enough samples for clustering, use single cluster
            labels = np.zeros(len(x), dtype=np.int8)
            labels[sample_mask] = 0
            return 1, labels

        # Test different numbers of clusters
        best_k = self.min_clusters
        best_score = np.inf
        best_gmm = None
        best_labels_valid = None

        for k in range(self.min_clusters, min(self.max_clusters + 1, n_samples // 3)):
            # Fit GMM with k components
            gmm = GaussianMixture(
                n_components=k,
                covariance_type='full',
                max_iter=100,
                n_init=10,
                init_params='kmeans',
                random_state=self.random_state
            )

            try:
                gmm.fit(X)
                labels_valid = gmm.predict(X)

                # Compute criterion
                if self.criterion == 'AIC':
                    score = gmm.aic(X)
                elif self.criterion == 'BIC':
                    score = gmm.bic(X)
                else:  # ICL
                    score = self._compute_icl(gmm, X, labels_valid)

                # Select model with lowest score
                if score < best_score:
                    best_score = score
                    best_k = k
                    best_gmm = gmm
                    best_labels_valid = labels_valid

            except Exception as e:
                # GMM fitting failed (e.g., singular covariance)
                # Skip this k value
                continue

        # If no valid model found, use single cluster
        if best_gmm is None:
            labels = np.zeros(len(x), dtype=np.int8)
            labels[sample_mask] = 0
            return 1, labels

        # Create full labels array (including excluded samples)
        labels = np.full(len(x), -9, dtype=np.int8)  # -9 = excluded
        labels[sample_mask] = best_labels_valid

        return best_k, labels

    def _compute_icl(self, gmm: GaussianMixture, X: np.ndarray,
                     labels: np.ndarray) -> float:
        """
        Compute Integrated Complete-data Likelihood (ICL) criterion.

        ICL = BIC + 2 * Entropy
        where Entropy = -ΣΣ z_ik * log(gamma_ik)

        Args:
            gmm: Fitted GaussianMixture model
            X: Data array (n_samples, 2)
            labels: Cluster assignments

        Returns:
            ICL score
        """
        # Get BIC
        bic = gmm.bic(X)

        # Compute entropy
        # Get posterior probabilities gamma_ik = P(k|x_i)
        gamma = gmm.predict_proba(X)

        # Hard assignments z_ik (0 or 1)
        n_samples = X.shape[0]
        k = gmm.n_components
        z = np.zeros((n_samples, k))
        for i in range(n_samples):
            z[i, labels[i]] = 1

        # Entropy = -ΣΣ z_ik * log(gamma_ik)
        # Add small constant to avoid log(0)
        entropy = -np.sum(z * np.log(gamma + 1e-10))

        # ICL = BIC + 2 * Entropy
        icl = bic + 2 * entropy

        return icl
