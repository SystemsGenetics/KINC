#include "PairWiseCluster.h"


/**
 * Constructor.
 */
PairWiseCluster::PairWiseCluster(PairWiseSet * pwset, char ** method, int num_methods) {

  // Initialize the class members.
  this->cluster_samples = NULL;
  this->cluster_size = 0;
  this->num_missing = 0;
  this->num_outliers = 0;
  this->num_threshold = 0;
  this->num_goutliers = 0;
  this->neighbor = NULL;
  this->pwsim = NULL;
  this->index = 0;
  this->pwset = pwset;
  this->method = method;
  this->num_methods = num_methods;

  this->pwsim = (PairWiseSimilarity **) malloc(sizeof(PairWiseSimilarity *) * num_methods);
}
/**
 * Desctructor.
 */
PairWiseCluster::~PairWiseCluster() {
  if (this->cluster_samples) {
    free(this->cluster_samples);
  }
  for (int i = 0; i < this->num_methods; i++) {
    if (this->pwsim[i]) {
      delete this->pwsim[i];
    }
  }
  delete this->pwsim;
}

void PairWiseCluster::doSimilarity(int min_obs) {

  // Create a new PairWiseSet for this cluster
  for (int i = 0; i < this->num_methods; i++) {
    if (strcmp(method[i], "sc") == 0) {
      SpearmanSimilarity * ssim = new SpearmanSimilarity(this->pwset, min_obs, this->cluster_samples);
      ssim->run();
      this->pwsim[i] = (PairWiseSimilarity *) ssim;
    }
    else if (strcmp(method[i], "pc") == 0) {
      PearsonSimilarity * psim = new PearsonSimilarity(this->pwset, min_obs, this->cluster_samples);
      psim->run();
      this->pwsim[i] = (PairWiseSimilarity *) psim;
    }
  //  else if (strcmp(method[i], "mi") == 0) {
  //    MISimilarity * msim = new MISimilarity(this->pwset, min_obs, this->cluster_samples);
  //    msim->run();
  //    this->pwsim[i] = (PairWiseSimilarity *) msim;
  //  }
  }
}
/**
 * Sets the PairWiseCluster samples.
 */
void PairWiseCluster::setClusterSamples(int * samples, bool from_clean) {

  this->cluster_samples = (int *) malloc(sizeof(int) * this->pwset->n_orig);

  // If the samples list was derived using the clean samples of the PWSet
  // then we need to stretch back out the samples to their original size
  // and mark missing samples with a 9.
  int k = 0;
  if (from_clean) {
    for (int i = 0; i < this->pwset->n_orig; i++) {
      if (this->pwset->samples[i] == 1) {
        this->cluster_samples[i] = samples[k];
        k++;
      }
      else {
        this->cluster_samples[i] = this->pwset->samples[i];
      }
      if (this->cluster_samples[i] == 1) {
        this->cluster_size++;
      }
      else if (this->cluster_samples[i] == 8) {
        this->num_outliers++;
      }
      else if (this->cluster_samples[i] == 6) {
        this->num_threshold++;
      }
      else if (this->cluster_samples[i] == 7) {
        this->num_goutliers++;
      }
      else if (this->cluster_samples[i] == 9) {
        this->num_missing++;
      }
    }
  }
  // The provided samples are not derived from the clean set so just copy them
  else {
    for (int i = 0; i < this->pwset->n_orig; i++) {
      this->cluster_samples[i] = samples[i];
      if (samples[i] == 1) {
        this->cluster_size++;
      }
    }
  }
}
/**
 *
 */
void PairWiseCluster::printCluster() {
  printf("%i\t%i\t", this->pwset->gene1, this->pwset->gene2);
  if (this->cluster_samples) {
    for (int i = 0; i < this->pwset->n_orig; i++) {
      printf("%i", this->cluster_samples[i]);
    }
  }
  printf("\n");
}


