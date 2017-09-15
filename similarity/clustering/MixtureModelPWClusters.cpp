#include "MixtureModelPWClusters.h"

/**
 * Constructor for MixMod class.
 *
 * @param EMatrix *ematrix;
 */
MixtureModelPWClusters::MixtureModelPWClusters(PairWiseSet *pwset, int min_obs,
    char ** method, int num_methods) {

  // Initialize some values.
  this->pwset = pwset;
  this->min_obs = min_obs;
  this->method = method;
  this->num_methods = num_methods;
  labels = NULL;
  pwcl = NULL;
  data = NULL;

  // Create the PairWiseClusterList object
  pwcl = new PairWiseClusterList(pwset);

  // Make sure we have the correct number of observations before preparing
  // for the comparision.
  if (this->pwset->n_clean < this->min_obs) {
    return;
  }

  // Create the Gaussian Data object and set the dataDescription object.
  this->data = (double **) malloc(sizeof(double *) * pwset->n_clean);
  for (int i = 0; i < pwset->n_clean ; i++) {
    this->data[i] = (double *) malloc(sizeof(double) * 2);
    this->data[i][0] = pwset->x_clean[i];
    this->data[i][1] = pwset->y_clean[i];
  }

}

/**
 * Destructor for the MixMod class.
 */
MixtureModelPWClusters::~MixtureModelPWClusters() {
  // Free the data array.
  if (data) {
    for (int i = 0; i < pwset->n_clean; i++) {
      free(data[i]);
    }
    free(data);
  }
  if (labels) {
    free(labels);
  }
  delete pwcl;
}

/**
 * Executes mixture model clustering.
 */
void MixtureModelPWClusters::run(char * criterion, int max_clusters) {
  int64_t lowest_cluster_num = 9;
  int found = 0;

  // Make sure we have the correct number of observations before performing
  // the comparision.
  if (this->pwset->n_clean < this->min_obs) {
    return;
  }

  // The data we are using is quantitative, so set the data type.
  XEM::GaussianData * gdata = new XEM::GaussianData(pwset->n_clean, 2, data);
  XEM::DataDescription dataDescription(gdata);

  // Set the number of clusters to be tested to 9.
  vector<int64_t> nbCluster;
  for (int i = 1; i <= max_clusters; i++) {
    nbCluster.push_back(i);
  }

  // Create the ClusteringInput object.
  XEM::ClusteringInput cInput(nbCluster, dataDescription);

  // Set the criterion
  cInput.removeCriterion(0);
  if (strcmp(criterion, "BIC") == 0) {
    cInput.addCriterion(XEM::BIC);
  }
  else if (strcmp(criterion, "ICL") == 0) {
    cInput.addCriterion(XEM::ICL);
  }
  else if (strcmp(criterion, "NEC") == 0) {
    cInput.addCriterion(XEM::NEC);
  }
  else if (strcmp(criterion, "CV") == 0) {
    cInput.addCriterion(XEM::CV);
  }
  else if (strcmp(criterion, "DCV") == 0) {
    cInput.addCriterion(XEM::DCV);
  }

  // Finalize input: run a series of sanity checks on it
  cInput.finalize();

  // Run XEM::ClusteringMain 5 times to find the iteration with the
  // lowest number of clusters.  We do this because the MixModLib can
  // detect different numbers of clusters depending on the random starting
  // point in the algorithm.
  for (int i = 0; i < 1; i++) {

    // Find clusters.
    XEM::ClusteringMain cMain(&cInput);
    cMain.run();

    // Create a new XEM::ClusteringOutput object
    XEM::ClusteringOutput * cOutput = cMain.getOutput();

    // Order the clusters using BIC.
    if (strcmp(criterion, "BIC") == 0) {
      cOutput->sort(XEM::BIC);
    }
    else if (strcmp(criterion, "ICL") == 0) {
      cOutput->sort(XEM::ICL);
    }
    else if (strcmp(criterion, "NEC") == 0) {
      cOutput->sort(XEM::NEC);
    }
    else if (strcmp(criterion, "CV") == 0) {
      cOutput->sort(XEM::CV);
    }
    else if (strcmp(criterion, "DCV") == 0) {
      cOutput->sort(XEM::DCV);
    }

    if (cOutput->atLeastOneEstimationNoError()) {
      found = 1;

      // Get the best XEM::ClusteringModelOutput
      XEM::ClusteringModelOutput* cMOutput = cOutput->getClusteringModelOutput().front();
      XEM::LabelDescription * ldescription = cMOutput->getLabelDescription();
      XEM::Label * label = ldescription->getLabel();
      int64_t * cLabels = label->getTabLabel();

      // Find the number of clusters
      int64_t num_clusters = 0;
      for (int i = 0; i < this->pwset->n_clean; i++) {
        if (cLabels[i] > num_clusters) {
          num_clusters = cLabels[i];
        }
      }

      // If this is the fewest number of clusters we've seen then copy
      // the labels over.
      if (num_clusters < lowest_cluster_num) {
        lowest_cluster_num = num_clusters;
        if (this->labels == NULL) {
          this->labels = (int64_t *) malloc(sizeof(int64_t) * this->pwset->n_clean);
        }
        for (int i = 0; i < this->pwset->n_clean; i++) {
          this->labels[i] = cLabels[i];
        }
      }
      delete cLabels;
    }
  }
  delete gdata;

  // If a cluster was not found then simply return.
  if (!found) {
    return;
  }

  // Iterate through the clusters in the set that was selected.  We are
  // done iterating through the clusters when we run out of labels in the
  // labels array.
  int cluster_num = 1;
  int cluster_samples[this->pwset->n_clean];
  bool done = false;
  while (!done) {
    // As a default set done to be true. If we find a label for the current
    // cluster number then it will get set to false and we can deal with the
    // cluster.
    done = true;

    // Prepare arrays for outlier detection and removal.
    float cx[this->pwset->n_clean];
    float cy[this->pwset->n_clean];
    for (int j = 0; j < this->pwset->n_clean; j++) {
      cx[j] = 0;
      cy[j] = 0;
    }

    // Build the samples array of 0's and 1's to indicate which
    // samples are in the cluster and which are not, also populate the
    // cx and cy arrays for outlier detection.
    int l = 0;
    for (int i = 0; i < this->pwset->n_clean; i++) {
      if (this->labels[i] == cluster_num) {
        cx[l] = this->pwset->x_clean[i];
        cy[l] = this->pwset->y_clean[i];
        l++;
        done = false;
        cluster_samples[i] = 1;
      }
      else {
        cluster_samples[i] = 0;
      }
    }

    // Discover any outliers for clusters with size >= min_obs
    Outliers * outliersCx = NULL;
    Outliers * outliersCy = NULL;
    outliersCx = outliers_iqr(cx, this->pwset->n_clean, 1.5);
    outliersCy = outliers_iqr(cy, this->pwset->n_clean, 1.5);

    // Remove any outliers
    for (int i = 0; i < this->pwset->n_clean; i++) {
      for (int j = 0; j < outliersCx->n; j++) {
        if (cluster_samples[i] == 1 && this->pwset->x_clean[i] == outliersCx->outliers[j]) {
          cluster_samples[i] = 8;
        }
      }
      for (int j = 0; j < outliersCy->n; j++) {
        if (cluster_samples[i] == 1 && this->pwset->y_clean[i] == outliersCy->outliers[j]) {
          cluster_samples[i] = 8;
        }
      }
    }
    if (outliersCx) {
      free(outliersCx->outliers);
      free(outliersCx);
    }
    if (outliersCy) {
      free(outliersCy->outliers);
      free(outliersCy);
    }

    // If we found samples with the current cluster_num then create a
    // cluster and add it to the list.
    if (!done) {
      PairWiseCluster * cluster = new PairWiseCluster(this->pwset, this->method, this->num_methods);
      cluster->setClusterSamples(cluster_samples, true);
      cluster->doSimilarity(this->min_obs);
//        cluster->printCluster();
      this->pwcl->addCluster(cluster);
    }
    cluster_num++;
  }
}
