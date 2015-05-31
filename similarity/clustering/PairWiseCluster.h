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
    PairWiseSimilarity * pwsim;
    // The number of samples that are missing from the cluster.
    int num_missing;

    // The next cluster in the linked list.
    PairWiseCluster * neighbor;
    // The index of this cluster in a PairWiseList
    int index;

  public:
    PairWiseCluster(PairWiseSet * pwset);
    ~PairWiseCluster();
    void setPWSimilarity(PairWiseSimilarity * pwsim);
    void setClusterSamples(int * samples, bool from_clean);
    void doSimilarity(const char * method, int min_obs);
    void printCluster();
};

#endif

