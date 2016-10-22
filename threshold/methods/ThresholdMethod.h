#ifndef _THRESHOLDMETHOD_
#define _THRESHOLDMETHOD_

#include <getopt.h>

#include "../../ematrix/EMatrix.h"

void print_threshold_usage();

/**
 * A base class for thresholding methods to inherit from.
 *
 * This class constructor will read arguments from the command-line, and
 * provides the getters and setters for commonly needed class members.
 */
class ThresholdMethod {

  protected:
    // The clustering method used.
    char * clustering;
    // The expression matrix object.
    EMatrix * ematrix;
    // The directory where the expression matrix is found
    char * bin_dir;
    // Specifies the correlation method that was used: pc, mi, sc
    char ** method;
    // Indicates the number of methods.
    int num_methods;
    // The method (e.g. sc, mi, pc) to use for thersholding
    char * th_method;
    // The index of the th_method in the methods array
    int th_method_index;


    // DATA FILTERS FOR CLUSTERED SIMILARITY DATA
    // ------------------------------------------
    // The maximum number of missing values in the comparision.
    int max_missing;
    // The minimum number of samples in a cluster.
    int min_cluster_size;
    // The maximum number of clusters a pair-wise comparision can have.
    int max_modes;
    // The mininum range of expression a cluster must have.
    int min_range;

    float ** parseScores(char * scores_str);


  public:
    ThresholdMethod(EMatrix *ematrix, char ** method, int num_methods,
        char * th_method, char * clustering,
        int min_cluster_size, int max_missing, int max_modes, int min_range);
    ~ThresholdMethod();

    // GETTERS
    // -------
    int getMaxMissing() { return max_missing; }
    int getMinClusterSize() { return min_cluster_size; }

    // TO BE IMPLEMENTED BY THE CHILD CLASS
    // ------------------------------------
    double findThreshold();
};

#endif
