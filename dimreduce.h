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
#include "error.h"

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
} PairWiseSet;

int do_dimreduce(int argc, char *argv[]);
void pairwise_reduce(double *a2, int x, double *b2, int y, int n2,
    EMatrix ematrix, CCMParameters params, int * kept, FILE * cf);
void write_reduced_ematrix_line(PairWiseSet pws, FILE * cf);
void print_dimreduce_usage();

#endif
