#ifndef _CLUSTERS_
#define _CLUSTERS_

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

/**
 * The PairWiseCluster struct contains the information about a cluster.
 *
 */
class PairWiseCluster {
  public:

    // An array of zeros and ones indicating which samples comprise the cluster
    int * samples;
    // The number of samples in the samples array.
    int num_samples;
    // The number of samples in the cluster
    int cluster_size;
    // The similarity score. Can be a Pearson's Correlation coefficient, or
    // a Spearman's coefficent, MI score, etc.
    double score;
    // The next cluster in the linked list.
    PairWiseCluster * neighbor;

    PairWiseCluster();
    ~PairWiseCluster() {};
};

/**
 * A set of clusters that belong to the same two genes.
 */
class PairWiseClusterSet {
  public:
    // The indecies of the genes in the input Expression matrix
    int gene1;
    int gene2;
    // The method used to perform the pair-wise comparision.
    char * method;
    // The number of clusters in the set.
    int num_clusters;

    // The head for the linked list of PairWiseClusters.
    PairWiseCluster * head;

    PairWiseClusterSet(int gene1, int gene2, const char * method);
    ~PairWiseClusterSet();

    void addCluster(PairWiseCluster * pwc);
};

/**
 *
 */
class PairWiseClusterWriter {
  private:
    // An array of file pointers.
    FILE ** fps;
    // Specifies the correlation method: pc, mi, sc
    char * method;
    // The prefix for the filename.
    char * fileprefix;
    // A unique id to differentiate between parallel executions.
    int id;

    // Opens and creates file pointers for all of the
    void openOutFiles();
    void closeOutFiles();
  public:
    // Constructor.
    PairWiseClusterWriter(const char * method, const char * fileprefix, int id);
    // Destructor.
    ~PairWiseClusterWriter();
    // Writes a PairWiseCluster to the proper file.
    void writeClusters(PairWiseClusterSet * pwcs);
};


#endif
