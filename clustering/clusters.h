#ifndef _CLUSTERS_
#define _CLUSTERS_

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

/**
 * The PairWiseClusters struct contains the information about a cluster. It is
 * also a list node (has a *next element) to allow these objects to be
 * strung together in a list.
 */
class SampleCluster {
  public:
    // An array of zeros and ones indicating which samples from the input file
    // are to be used for the comparison.
    int * samples;
    // The number of samples in the samples array.
    int num_samples;
    // The index of the first gene used in the comparison.
    int gene1;
    // The index of the second gene used in the comparison.
    int gene2;
    // The number of samples in the cluster
    int cluster_size;
    // The Pearson's Correlation Coefficient for the cluster
    double pcc;

    SampleCluster * next;

    SampleCluster();
    ~SampleCluster();

    void addSampleCluster(SampleCluster * newc);

};

class SampleClusterWriter {
  private:
    // An array of file pointers.
    FILE ** fps;

    void openOutFiles(char * method, char * fileprefix, int mpi_id);
    void closeOutFiles();
  public:
    SampleClusterWriter();
    ~SampleClusterWriter();
    void writeSampleCluster(SampleCluster pwc);
};


//FILE ** open_output_files(CCMParameters params, int mpi_id);
//void close_output_files(FILE** fps);
//void write_pairwise_cluster_samples(PairWiseClusters * pwc, FILE ** fps);
//// Functions for working with the PairWiseClusters list
//void free_pairwise_cluster_list(PairWiseClusters * head);
//PairWiseClusters * new_pairwise_cluster_list();
//void add_pairwise_cluster_list(PairWiseClusters **head, PairWiseClusters *newc);
//void update_pairwise_cluster_samples(int * parent_samples, int n, PairWiseClusters * head);


#endif
