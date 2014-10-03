#ifndef _PREPROCESS_
#define _PREPROCESS_

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <time.h>
#include <string.h>
#include "similarity.h"
#include "stats/royston.h"

int do_preprocess(int argc, char *argv[]);
void print_preprocess_usage();

#endif
