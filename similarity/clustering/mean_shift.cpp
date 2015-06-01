#include "mean_shift.h"

/**
 * @param double * a2
 *   A row from the ematrix with missing values removed
 * @param int x
 *   The index of a2 in the ematrix
 * @param double * b2
 *   A row from the ematrix with missing values removed.
 * @param in y
 *   The index of b2 in the ematrix
 * @param int n2
 *   The size of both a2 and b2.
 * @param Ematrix ematrix
 *   The two-dimensional ematrix where rows are genes and columns are samples.
 * @param CCMParameters params
 *   The parameters provided to the program
 * @param float bw
 *   The bandwith argument for mean shift clustering.  Default should be 0.075.
 * @param int level
 *   An integer indicating the recursion level. It should always be set to
 *   zero by the caller of this function.
 */
PairWiseClusters * clustering(double *a2, int x, double *b2, int y, int n2,
    EMatrix * ematrix, CCMParameters params, float bw, int level) {

  // Variables used for looping.
  int k, j;
  int nkept;
  int * ckept;

  PairWiseClusters * result = new_pairwise_cluster_list();

  // Perform bandwidth selection
  double * selected = meanshift_coverage2D(a2, b2, n2);

  // Perform mean shift clustering (MSC)
  MeanShiftClusters * clusters;
  clusters = meanshift2D(a2, b2, n2, selected[0]);

  // Count the number of clusters that are larger than min_obs
  int num_large = 0;
  for(k = 0; k < clusters->num_clusters; k++) {
    if (clusters->sizes[k] >= params.min_obs) {
      num_large++;
    }
  }

  // Iterate through all of the clusters.
  for(k = 0; k < clusters->num_clusters; k++) {

    // For easier access create variables for the current cluster size.
    int size = clusters->sizes[k];

    // Create separate vectors for the x and y coordinates.
    double *cx = (double *) malloc(sizeof(double) * size);
    double *cy = (double *) malloc(sizeof(double) * size);
    for (j = 0; j < size; j++) {
      cx[j] = 0;
      cy[j] = 0;
    }
    // Add the values to the cx & cy arrays
    int l = 0;
    for (j = 0; j < n2; j++) {
      if (clusters->cluster_label[j] == k + 1) {
        cx[l] = a2[j];
        cy[l] = b2[j];
        l++;
      }
    }

    // Discover any outliers for clusters with size >= min_obs
    Outliers * outliersCx = NULL;
    Outliers * outliersCy = NULL;
    if (clusters->sizes[k] >= params.min_obs) {
      outliersCx = outliers_iqr(cx, size, 1.5);
      outliersCy = outliers_iqr(cy, size, 1.5);
    }

    // Create an array of the kept samples for this cluster. Don't include
    // any samples whose coordinates are considered outliers or who are not
    // in the cluster.
    ckept = (int *) malloc(sizeof(int) * n2);
    for (l = 0; l < n2; l++) {
      ckept[l] = 0;
    }
    nkept = 0;
    for (l = 0; l < n2; l++) {
      // Is this sample is in the cluster? if so then also make sure it's
      // not an outlier.
      if (clusters->cluster_label[l] == k + 1) {
        // Iterate through the outlier points and compare to this
        // sapmle's points.  If there is a match then mark as an outlier
        // and exclude it from the samples that are kept.  First check the
        // x coordinate
        int is_outlier = 0;
        if (outliersCx) {
          for (j = 0; j < outliersCx->n; j++) {
            if (a2[l] == outliersCx->outliers[j]) {
              is_outlier = 1;
              break;
            }
          }
        }
        // Second check the y coordinate.
        if (outliersCy) {
          for (j = 0; j < outliersCy->n; j++) {
            if (b2[l] == outliersCy->outliers[j]) {
              is_outlier = 1;
              break;
            }
          }
        }
        // if it's not an outlier then keep it.
        if (!is_outlier) {
          ckept[l] = 1;
          cx[nkept] = a2[l];
          cy[nkept] = b2[l];
          nkept++;
        }
      }
    } // end for (l = 0; l < n2; l++) ...

    // Free the outlier objects.
    if (outliersCx) {
      free(outliersCx->outliers);
      free(outliersCx);
    }
    if (outliersCy) {
      free(outliersCy->outliers);
      free(outliersCy);
    }

    // Now after we have removed outliers and non-cluster samples,
    // makes sure we still have the minimum observations.
    if (nkept >= params.min_obs) {

      // If the recursion level is greater than zero the we've already clustered
      // at least once.  When we reach this point we've clustered again. If
      // we only have a single cluster at this point then we can't cluster
      // anymore and we can perform correlation.
      if (level > 0 && num_large == 1) {
        float corr = NAN;
        if (strcmp(params.method, "sc") == 0) {
          // Initialize the workspace needed for Spearman's calculation.
          double workspace[2 * params.rows];
          corr = gsl_stats_spearman(cx, 1, cy, 1, nkept, workspace);
        }
        if (strcmp(params.method, "pc") == 0) {
          corr = gsl_stats_correlation(cx, 1, cy, 1, nkept);
        }
        PairWiseClusters * newc = new_pairwise_cluster_list();
        newc->gene1 = x;
        newc->gene2 = y;
        newc->num_samples = n2;
        newc->samples = ckept;
        newc->next = NULL;
        newc->cluster_size = nkept;
        newc->pcc = corr;
        add_pairwise_cluster_list(&result, newc);
      }

      // If this is the first level of clustering then we want to cluster
      // again, but only if there is one or more large cluster.
      // We do this because we haven't yet settled on a distinct
      // "expression mode".
      else {
        PairWiseClusters * children = clustering(cx, x, cy, y, nkept, ematrix, params, 0.09, level + 1);
        update_pairwise_cluster_samples(ckept, n2, children);
        add_pairwise_cluster_list(&result, children);
        // here we need to free ckept because it never got used in a
        // PairWiseClusters instance
        free(ckept);
      }
    }
    // If the cluster is too small then just add it without
    // correlation analysis or further sub clustering.
    else {
      PairWiseClusters * newc = new_pairwise_cluster_list();
      newc->gene1 = x;
      newc->gene2 = y;
      newc->num_samples = n2;
      newc->samples = ckept;
      newc->next = NULL;
      newc->cluster_size = nkept;
      newc->pcc = NAN;
      add_pairwise_cluster_list(&result, newc);
    }
    free(cx);
    free(cy);
  }

  // free the memory
  free_msc(clusters);


  return result;
}


