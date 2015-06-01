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
class MixtureModelClustering : public PairWiseClustering {
  private:
    // The criterion model. E.g.  BIC, ICL, NEC, CV, DCV.
    char * criterion;
    // The maximum number of clusters to allow per comparision.
    int max_clusters;
    // The similarity method
    char * method;

  public:

    MixtureModelClustering(EMatrix *ematrix, int min_obs, int num_jobs,
        int job_index, char * method, char * criterion, int max_clusters);
    ~MixtureModelClustering();

    // Returns usage help instructions.
    void getUsage();
    // Performs pair-wise mixture model clustering of the entire ematrix.
    void run();

};

#endif
