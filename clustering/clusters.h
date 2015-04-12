#ifndef _CLUSTERS_
#define _CLUSTERS_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <libgen.h>
#include <string.h>
#include <math.h>
#include "../similarity.h"
#include "../similarity/spearman.h"
#include "../similarity/pearson.h"
#include "../similarity/bspline_mi.h"


/**
 * The PairWiseCluster class contains the information about a cluster.
 *
 */
class PairWiseCluster {
  friend class PairWiseClusterList;
  friend class PairWiseClusterWriter;

  public:
    // An array of zeros and ones indicating which samples comprise the cluster.
    // this array will be the same size as the n_orig of the PairWiseSet.
    int * cluster_samples;
    // The number of samples in the cluster
    int cluster_size;
    // The next cluster in the linked list.
    PairWiseCluster * neighbor;
    // The object containing the pair of genes/probesets to perform clustering.
    PairWiseSet * pwset;
    // A similarity function object.
    PairWiseSimilarity * pwsim;

  public:
    PairWiseCluster(PairWiseSet * pwset);
    ~PairWiseCluster();
    void setPWSimilarity(PairWiseSimilarity * pwsim);
    void setClusterSamples(int * samples, bool from_clean);
    void doSimilarity(const char * method, int min_obs);
    void printCluster();
};

/**
 * A set of clusters that belong to the same two genes.
 */
class PairWiseClusterList {
  friend class PairWiseClusterWriter;
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
    ~PairWiseClusterList() {};

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
    PairWiseClusterWriter(char * method, char * fileprefix, int id);
    // Destructor.
    ~PairWiseClusterWriter();
    // Writes a PairWiseCluster to the proper file.
    void writeClusters(PairWiseClusterList *pwcl, int gene1, int gene2);
};


#endif
