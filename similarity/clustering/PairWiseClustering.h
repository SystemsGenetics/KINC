#ifndef _CLUSTERING_
#define _CLUSTERING_

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <libgen.h>
#include <signal.h>

#include "../../ematrix/EMatrix.h"
#include "../../general/error.h"
#include "../../general/misc.h"
#include "PairWiseCluster.h"

/**
 * A base class for pair-wise clustering of samples.
 */
class PairWiseClustering {

  protected:
    // The expression matrix object.
    EMatrix * ematrix;
    // The minimum number of observations to calculate correlation.
    int min_obs;
    // The total number of jobs that will be run at once.
    int num_jobs;
    // The index of this job within the total jobs.  Must be
    // between 1 and num_jobs (no zero index).
    int job_index;

  public:
    PairWiseClustering(EMatrix *ematrix, int min_obs, int num_jobs, int job_index);
    ~PairWiseClustering();

    // Returns the minimum number of observations required to form a cluster.
    int getMinObs() { return min_obs; }
    // Returns the number of clustering jobs being run.
    int getNumJobs() { return num_jobs; }
    // Retuns the job index for this thread.
    int getJobIndex() { return job_index; }
    // Returns the usage instructions.
    char * getUsage();

    // Perform pair-wise mixture model clustering for every gene pair.
    // To be implemented by class children.
    void run() {}

};


#endif
