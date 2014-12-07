#ifndef _DIMREDUCE_
#define _DIMREDUCE_

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include "similarity.h"
#include "stats/royston.h"
#include "stats/meanshift.h"
#include "stats/outlier.h"
#include "error.h"

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
  // The cluster label.  Set to 0 if no clustering was performed.
  int cluster_label;

  struct PairWiseClusters * next;
} PairWiseClusters;


// Primary function for this file
int do_dimreduce(int argc, char *argv[]);
void print_dimreduce_usage();

// Functions for working with the PairWiseClusters list
PairWiseClusters * new_pairiwse_cluster_list();
void add_pairiwse_cluster_list(PairWiseClusters *head, PairWiseClusters *new);
void update_pairwise_cluster_samples(int * ckept, int nkept, PairWiseClusters *new);
void write_pairwise_cluster_samples(PairWiseClusters *pws, FILE * cf);

// Function to peform the clustering
PairWiseClusters clustering(double *a2, int x, double *b2, int y, int n2,
    EMatrix ematrix, CCMParameters params);


#endif
