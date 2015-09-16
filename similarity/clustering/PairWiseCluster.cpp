#include "PairWiseCluster.h"


/**
 * Constructor.
 */
PairWiseCluster::PairWiseCluster(PairWiseSet * pwset) {

  // Initialize the class members.
  this->cluster_samples = NULL;
  this->cluster_size = 0;
  this->num_missing = 0;
  this->num_outliers = 0;
  this->neighbor = NULL;
  this->pwsim = NULL;
  this->index = 0;
  this->pwset = pwset;
}
/**
 * Desctructor.
 */
PairWiseCluster::~PairWiseCluster() {
  free(this->cluster_samples);
  if (this->pwsim) {
    delete this->pwsim;
  }
}

void PairWiseCluster::doSimilarity(const char * method, int min_obs) {
  // Create a new PairWiseSet for this cluster
  if (strcmp(method, "sc") == 0) {
    SpearmanSimilarity * ssim = new SpearmanSimilarity(this->pwset, min_obs, this->cluster_samples);
    ssim->run();
    this->pwsim = (PairWiseSimilarity *) ssim;
  }
  else if (strcmp(method, "pc") == 0) {
    PearsonSimilarity * psim = new PearsonSimilarity(this->pwset, min_obs, this->cluster_samples);
    psim->run();
    this->pwsim = (PairWiseSimilarity *) psim;
  }
  else if (strcmp(method, "mi") == 0) {
//    MISimilarity * msim = new MISimilarity(this->pwset, min_obs, this->cluster_samples);
//    msim->run();
//    this->pwsim = (PairWiseSimilarity *) msim;
  }
}
/**
 * Sets the PairWiseCluster samples.
 */
void PairWiseCluster::setClusterSamples(int * samples, bool from_clean) {

  this->cluster_samples = (int *) malloc(sizeof(int) * this->pwset->n_orig);

  // If the samples list was derived using the clean samples of the PWSet
  // then we need strech back out the samples to their original size
  // and marking missing samples with a 9.
  int k = 0;
  if (from_clean) {
    for (int i = 0; i < this->pwset->n_orig; i++) {
      if (this->pwset->samples[i] == 1) {
        this->cluster_samples[i] = samples[k];
        if (samples[k] == 1) {
          this->cluster_size++;
        }
        if (samples[k] == 8) {
          this->num_outliers++;
        }
        k++;
      }
      else {
        this->cluster_samples[i] = 9;
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


