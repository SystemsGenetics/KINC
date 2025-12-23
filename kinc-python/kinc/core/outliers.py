"""
Outlier detection using Tukey's method (IQR-based).
"""

import numpy as np
from typing import Tuple


def remove_outliers(x: np.ndarray, y: np.ndarray,
                    labels: np.ndarray, cluster: int) -> np.ndarray:
    """
    Remove outliers from a specific cluster using Tukey's method.

    Uses the interquartile range (IQR) method:
    - Outliers are values outside [Q1 - 1.5*IQR, Q3 + 1.5*IQR]

    Args:
        x: Expression values for gene X
        y: Expression values for gene Y
        labels: Cluster labels for each sample
        cluster: Which cluster to process

    Returns:
        Updated labels with outliers marked as -7 or -8
    """
    # Get samples in this cluster
    mask = labels == cluster
    x_cluster = x[mask]
    y_cluster = y[mask]

    if len(x_cluster) < 4:  # Need at least 4 points for quartiles
        return labels

    # Compute quartiles and IQR for x
    q1_x, q3_x = np.percentile(x_cluster, [25, 75])
    iqr_x = q3_x - q1_x
    lower_x = q1_x - 1.5 * iqr_x
    upper_x = q3_x + 1.5 * iqr_x

    # Compute quartiles and IQR for y
    q1_y, q3_y = np.percentile(y_cluster, [25, 75])
    iqr_y = q3_y - q1_y
    lower_y = q1_y - 1.5 * iqr_y
    upper_y = q3_y + 1.5 * iqr_y

    # Identify outliers (outside bounds in either dimension)
    outliers_x = (x_cluster < lower_x) | (x_cluster > upper_x)
    outliers_y = (y_cluster < lower_y) | (y_cluster > upper_y)
    outliers = outliers_x | outliers_y

    # Update labels (mark outliers with special value)
    # Make a copy to avoid modifying original
    new_labels = labels.copy()
    cluster_indices = np.where(mask)[0]
    outlier_indices = cluster_indices[outliers]

    # Determine outlier code based on current label value
    # -7 for pre-clustering outliers, -8 for post-clustering outliers
    outlier_code = -7 if cluster == 0 else -8
    new_labels[outlier_indices] = outlier_code

    return new_labels


def remove_outliers_all_clusters(x: np.ndarray, y: np.ndarray,
                                 labels: np.ndarray,
                                 n_clusters: int) -> np.ndarray:
    """
    Remove outliers from all clusters.

    Args:
        x: Expression values for gene X
        y: Expression values for gene Y
        labels: Cluster labels for each sample
        n_clusters: Number of clusters

    Returns:
        Updated labels with outliers marked
    """
    new_labels = labels.copy()

    for k in range(n_clusters):
        new_labels = remove_outliers(x, y, new_labels, cluster=k)

    return new_labels
