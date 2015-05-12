#ifndef _DIMREDUCE_
#define _DIMREDUCE_

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <libgen.h>
#include <signal.h>

#include "ematrix.h"
//#include "stats/royston.h"
#include "stats/meanshift.h"
//#include "stats/outlier.h"
#include "error.h"
#include "misc.h"
#include "clustering/mixmod.h"
#include "similarity.h"
#include "clustering/clusters.h"

/**
 * This class holds the arguments used as input to the dimreduce program.
 */
class DRArgs {

  private:
    // Indicates if missing values should be ignored in the EMatrix file.
    int omit_na;
    // Indicates if headers are present in the input EMatrix file.
    int headers;
    // The number of rows in the expression matrix.
    int rows;
    // The number of columns in the expression matrix.
    int cols;
    // The input file name
    char *infilename;
    // Specifies the value that represents a missing value.
    char *na_val;
    // Specifies the transformation function: log2, none.
    char func[10];
    // Specifies the correlation method: pc, mi, sc
    char method[10];
    // The minimum number of observations to calculate correlation.
    int min_obs;
    // Set to 1 to perform log10 transformation.
    int do_log10;
    // Set to 1 to perform log2 transformation.
    int do_log2;
    // Set to 1 to perform log transformation.
    int do_log;
    // The input filename without the prefix.
    char * fileprefix;
    // The total number of jobs that will be run at once.
    int num_jobs;
    // The index of this job within the total jobs.  Must be
    // between 1 and num_jobs (no zero index)
    int job_index;

    // Variables for mean shift clustering
    double msc_bw1;
    double msc_bw2;

    // Variables for mixture models
    // The criterion model. E.g.  BIC, ICL, NEC, CV, DCV
    char criterion[4];
    // Max clusters
    int max_clusters;

  public:
    DRArgs(int argc, char *argv[]);
    ~DRArgs();

    // Getters
    int getOmitNA() { return omit_na; }
    int getHasHeaders() { return headers; }
    int getNumRows() { return rows; }
    int getNumCols() { return cols; }
    char * getInfileName() { return infilename; }
    char * getNAval() { return na_val; }
    char * getTransformFunction() { return func; }
    char * getCorMethod() { return method; }
    int getMinObs() { return min_obs; }
    int getDoLog10() { return do_log10; }
    int getDoLog2() { return do_log2; }
    int getDoLog() { return do_log; }
    int getNumJobs() { return num_jobs; }
    int getJobIndex() { return job_index; }
    char * getFilePrefix() { return fileprefix; }

    char * getCriterion() { return criterion; }
    int getMaxClusters() { return max_clusters; }
};

// Primary function for this file
//int do_dimreduce(int argc, char *argv[], int mpi_id, int mpi_num_procs);
int do_dimreduce(int argc, char *argv[]);
void print_dimreduce_usage();


#endif
