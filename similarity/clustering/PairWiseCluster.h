#ifndef _PAIRWISECLUSTER_
#define _PAIRWISECLUSTER_

#include "../PairWiseSet.h"
#include "../PairWiseSimilarity.h"
#include "../SpearmanSimilarity.h"
#include "../PearsonSimilarity.h"
#include "../MISimilarity.h"

/**
 * The PairWiseCluster class contains the information about a cluster.
 *
 */
class PairWiseCluster {

  public:
    // An array of zeros and ones indicating which samples comprise the cluster.
    // this array will be the same size as the n_orig of the PairWiseSet.
    int * cluster_samples;
    // The number of samples in the cluster
    int cluster_size;
    // The object containing the pair of genes/probesets to perform clustering.
    PairWiseSet * pwset;
    // A similarity function object.
    PairWiseSimilarity ** pwsim;
    // The number of samples that are missing from the cluster.
    int num_missing;
    // The number of samples removed as outliers in the cluster.
    int num_outliers;
    // The number of samples removed as outliers from the set before clustering.
    int num_goutliers;
    // The number of samples excluded because values were below expression threshold.
    int num_threshold;
    // The pair-wise similarity methods to be performed.
    char ** method;
    // The number of similarity methods to perform.
    int num_methods;

    // The next cluster in the linked list.
    PairWiseCluster * neighbor;
    // The index of this cluster in a PairWiseList
    int index;

  public:
    PairWiseCluster(PairWiseSet * pwset, char ** method, int num_methods);
    ~PairWiseCluster();
    void setPWSimilarity(PairWiseSimilarity * pwsim);
    void setClusterSamples(int * samples, bool from_clean);
    void doSimilarity(int min_obs);
    void printCluster();

    PairWiseSimilarity * getPWSimilarity(int index) {return this->pwsim[index]; }
};

#endif

