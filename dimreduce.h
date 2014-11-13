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
  int * samples;
  int num_samples;
  char * gene1;
  char * gene2;
} PairWiseSet;

int do_dimreduce(int argc, char *argv[]);
void find_valid_comps(double *a2, double *b2, int n2, EMatrix ematrix, CCMParameters params);
void print_dimreduce_usage();

#endif
