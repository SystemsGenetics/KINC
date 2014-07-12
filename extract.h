#ifndef _EXTRACT_
#define _EXTRACT_

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <limits.h>
#include <dirent.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <libgen.h>


typedef struct{

  char * infilename;          // The input expression matrix file name.
  char * inputDir;

  char gene_names_file[1024];
  char fileprefix[1024];      // The input filename without the prefix.
  char method[10];            // Specifies the method: cor, mi.
  int x_coord;
  int y_coord;
  float th;                   // The threshold for creating the network.

  // the number of genes, n, in the n x n similarity matrix
  int numGenes;

  // the number of lines per bin file
  int numLinesPerFile;

  // holds file handles for all of the binary files.
  FILE * files[50];
  int num_files; // indicates the number of files in the files array

} NetParameters;

// Prototypes

int do_extract(int argc, char *argv[]);
void open_bin_files(NetParameters params);
void close_bin_files(NetParameters params);
void print_extract_usage();

#endif
