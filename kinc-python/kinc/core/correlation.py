"""
Correlation functions (Pearson and Spearman).
"""

import numpy as np
from typing import Optional


def pearson_correlation(x: np.ndarray, y: np.ndarray,
                       labels: Optional[np.ndarray] = None,
                       cluster: int = 0,
                       min_samples: int = 30) -> float:
    """
    Compute Pearson correlation coefficient for a gene pair.

    Args:
        x: Expression values for gene X
        y: Expression values for gene Y
        labels: Cluster labels (None = all samples in one cluster)
        cluster: Which cluster to compute correlation for
        min_samples: Minimum number of samples required

    Returns:
        Pearson correlation coefficient or NaN if insufficient samples
    """
    # Select samples for this cluster
    if labels is not None:
        mask = labels == cluster
        x_cluster = x[mask]
        y_cluster = y[mask]
    else:
        x_cluster = x
        y_cluster = y

    # Check minimum sample requirement
    n = len(x_cluster)
    if n < min_samples:
        return np.nan

    # Remove NaN values
    valid = ~(np.isnan(x_cluster) | np.isnan(y_cluster))
    x_valid = x_cluster[valid]
    y_valid = y_cluster[valid]

    if len(x_valid) < min_samples:
        return np.nan

    # Compute Pearson correlation using efficient formula
    # r = (n*Σxy - Σx*Σy) / sqrt((n*Σx² - (Σx)²) * (n*Σy² - (Σy)²))
    n = len(x_valid)
    sum_x = np.sum(x_valid)
    sum_y = np.sum(y_valid)
    sum_xx = np.sum(x_valid * x_valid)
    sum_yy = np.sum(y_valid * y_valid)
    sum_xy = np.sum(x_valid * y_valid)

    numerator = n * sum_xy - sum_x * sum_y
    denominator_x = n * sum_xx - sum_x * sum_x
    denominator_y = n * sum_yy - sum_y * sum_y

    if denominator_x <= 0 or denominator_y <= 0:
        return np.nan

    denominator = np.sqrt(denominator_x * denominator_y)
    if denominator == 0:
        return np.nan

    correlation = numerator / denominator

    # Clamp to [-1, 1] to handle floating point errors
    correlation = np.clip(correlation, -1.0, 1.0)

    return correlation


def spearman_correlation(x: np.ndarray, y: np.ndarray,
                        labels: Optional[np.ndarray] = None,
                        cluster: int = 0,
                        min_samples: int = 30) -> float:
    """
    Compute Spearman rank correlation coefficient for a gene pair.

    Args:
        x: Expression values for gene X
        y: Expression values for gene Y
        labels: Cluster labels (None = all samples in one cluster)
        cluster: Which cluster to compute correlation for
        min_samples: Minimum number of samples required

    Returns:
        Spearman correlation coefficient or NaN if insufficient samples
    """
    # Select samples for this cluster
    if labels is not None:
        mask = labels == cluster
        x_cluster = x[mask]
        y_cluster = y[mask]
    else:
        x_cluster = x
        y_cluster = y

    # Check minimum sample requirement
    if len(x_cluster) < min_samples:
        return np.nan

    # Remove NaN values
    valid = ~(np.isnan(x_cluster) | np.isnan(y_cluster))
    x_valid = x_cluster[valid]
    y_valid = y_cluster[valid]

    if len(x_valid) < min_samples:
        return np.nan

    # Compute ranks (scipy's rankdata with average method for ties)
    from scipy.stats import rankdata
    x_ranks = rankdata(x_valid, method='average')
    y_ranks = rankdata(y_valid, method='average')

    # Apply Pearson correlation to ranks
    return pearson_correlation(x_ranks, y_ranks, min_samples=min_samples)


def compute_correlations(x: np.ndarray, y: np.ndarray,
                         labels: np.ndarray,
                         n_clusters: int,
                         method: str = 'pearson',
                         min_samples: int = 30) -> np.ndarray:
    """
    Compute correlations for all clusters in a gene pair.

    Args:
        x: Expression values for gene X
        y: Expression values for gene Y
        labels: Cluster labels for each sample
        n_clusters: Number of clusters
        method: 'pearson' or 'spearman'
        min_samples: Minimum samples per cluster

    Returns:
        Array of correlations for each cluster
    """
    correlations = np.zeros(n_clusters, dtype=np.float32)

    corr_func = pearson_correlation if method == 'pearson' else spearman_correlation

    for k in range(n_clusters):
        correlations[k] = corr_func(x, y, labels, cluster=k, min_samples=min_samples)

    return correlations
