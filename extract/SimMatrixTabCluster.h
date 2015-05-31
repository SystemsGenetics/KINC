#ifndef _SIMMATRIXTABCLUSTER_
#define _SIMMATRIXTABCLUSTER_

#include <sys/stat.h>
#include "SimilarityMatrix.h"

/**
 * Class for extract a network from a tab-delimited cluster file.
 */
class SimMatrixTabCluster : public SimilarityMatrix {

  private:
    // The number of jobs that were used to construct the files.
    int num_jobs;


    // Discovers the number of jobs used to generate the clustering files.
    void getNumJobs();
  public:
    // Constructur.
    SimMatrixTabCluster(int argc, char *argv[]);
    // Destructor.
    ~SimMatrixTabCluster();
    // Retrieves the set of edges that match the given filtering parameters.
    // The user must have provided a threshold value.
    void writeNetwork();
    // Retrieves the similarity value for the given filtering paramters.
    // The user must have provided an x and y coordiante.
    void getPosition();
};

#endif
