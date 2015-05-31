#ifndef _MIXMODPWCLUSTERING_
#define _MIXMODPWCLUSTERING_

#include "../SpearmanSimilarity.h"
#include "PairWiseCluster.h"
#include "PairWiseClusterList.h"
#include "PairWiseClusterWriter.h"
#include "MixtureModelPWClusters.h"

/**
 * A class for performing mixuture models for an entire ematrix.
 */
class MixtureModelPWClustering : public PairWiseClustering {
  private:
    // The criterion model. E.g.  BIC, ICL, NEC, CV, DCV.
    char criterion[4];
    // The maximum number of clusters to allow per comparision.
    int max_clusters;

  public:

    MixtureModelPWClustering(int argc, char *argv[], EMatrix * ematrix);
    ~MixtureModelPWClustering();

    // Returns usage help instructions.
    void getUsage();
    // Performs pair-wise mixture model clustering of the entire ematrix.
    void run();

};

#endif
