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

int do_dimreduce(int argc, char *argv[]);
void print_dimreduce_usage();

#endif
