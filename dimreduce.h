#ifndef _DIMREDUCE_
#define _DIMREDUCE_

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include <mpi.h>
#include "similarity.h"
#include "stats/royston.h"
#include "stats/meanshift.h"
#include "stats/outlier.h"
#include "error.h"
#include "misc.h"

/**
 * The PairWiseClusters struct contains the information about a cluster. It is
 * also a list node (has a *next element) to allow these objects to be
 * strung together in a list.
 */
 typedef struct {
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
  // The Pearsons Correlation Coeeficient for the cluster
  double pcc;

  struct PairWiseClusters * next;
} PairWiseClusters;


// Primary function for this file
int do_dimreduce(int argc, char *argv[], int mpi_id, int mpi_num_procs);
void print_dimreduce_usage();

// Functions for working with the PairWiseClusters list
void free_pairwise_cluster_list(PairWiseClusters * head);
PairWiseClusters * new_pairwise_cluster_list();
void add_pairwise_cluster_list(PairWiseClusters *head, PairWiseClusters *new);
void write_pairwise_cluster_samples(PairWiseClusters * pwc, FILE ** fps);
void update_pairwise_cluster_samples(int * parent_samples, int n, PairWiseClusters * head);
FILE ** open_output_files(CCMParameters params, int mpi_id);
void close_output_files(FILE** fps);

// Function to peform the clustering
PairWiseClusters * clustering(double *a2, int x, double *b2, int y, int n2,
    EMatrix ematrix, CCMParameters params, float bw, int level);


#endif
