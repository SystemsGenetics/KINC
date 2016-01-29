#ifndef _MIXMODPWCLUSTERING_
#define _MIXMODPWCLUSTERING_

#include "PairWiseCluster.h"
#include "PairWiseClusterList.h"
#include "PairWiseClusterWriter.h"
#include "MixtureModelPWClusters.h"

// This structure is used for containing the set of genes for filtering.
struct geneFilter {
  char * filename;
  int * indicies;
  int num_genes;
};

/**
 * A class for performing mixuture models for an entire ematrix.
 */
class MixtureModelClustering : public PairWiseClustering {
  private:
    // The criterion model. E.g.  BIC, ICL, NEC, CV, DCV.
    char * criterion;
    // The maximum number of clusters to allow per comparision.
    int max_clusters;
    // The similarity method.
    char ** method;
    // Indicates the number of methods.
    int num_methods;
    // The threshold for expression values.
    double threshold;
    // Filter sets.
    geneFilter *set1, *set2;
    // The minimum similarity to write out.
    double * min_sim;

  public:

    MixtureModelClustering(EMatrix *ematrix, int min_obs, int num_jobs,
        int job_index, char ** method, int num_methods, char * criterion,
        int max_clusters, double threshold, geneFilter * set1,
        geneFilter * set2, double * min_sim);
    ~MixtureModelClustering();

    // Returns usage help instructions.
    void getUsage();
    // Performs pair-wise mixture model clustering of the entire ematrix.
    void run();

};

#endif
