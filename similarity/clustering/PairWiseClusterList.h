#ifndef _PWCLUSTERLIST_
#define _PWCLUSTERLIST_

#include "PairWiseCluster.h"

/**
 * A set of clusters that belong to the same two genes.
 */
class PairWiseClusterList {
  private:
    // The object containing the pair of genes/probesets to perform clustering.
    PairWiseSet * pwset;
    // When new clusters are added to the set, their samples may or may not
    // be the same size as the samples in this cluster set.  This is
    // because clustering may occur on subsets of samples as is the case when
    // samples are removed because of missing values. This function will
    // expand the samples array of the PairWiseCluster object to be the
    // same size as the set.
    void updateClusterSamples(PairWiseCluster * pwc);
  public:
    // The number of clusters in the set.
    int num_clusters;

    // The head for the linked list of PairWiseClusters.
    PairWiseCluster * head;

    PairWiseClusterList(PairWiseSet * pwset);
    ~PairWiseClusterList();

    void addCluster(PairWiseCluster * pwc);
};

#endif
