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
  char * genes_file;

  char fileprefix[1024]; // The input filename without the prefix.
  char method[10];       // Specifies the method: cor, mi.
  int x_coord;           // The user-specified x coordinate to retrieve
  int y_coord;           // The user-specified y coordinate to retrieve
  char * gene1;          // The user-specified gene1 of pair-wise comparision to retrieve
  char * gene2;          // the user-specified gene2 of pair-wise comparision to retrieve
  float th;              // The threshold for creating the network.

  // holds file handles for all of the binary files.
  FILE * files[50];  // an array of file pointers to the bin files
  int numLines[50];  // the number of lines per bin file
  int numGenes;      // the number of genes in the similarity matrix
  int num_files;     // indicates the number of files in the files array

  char ** genes;     // an array containing the names of genes, in order, for the smatrix

} NetParameters;

// Prototypes

int do_extract(int argc, char *argv[]);
void get_edges(NetParameters params);
void get_position(NetParameters params);
void open_bin_files(NetParameters *params);
void close_bin_files(NetParameters params);
void get_gene_names(NetParameters *params);
void get_gene_coords(NetParameters *params);
void print_extract_usage();

#endif
