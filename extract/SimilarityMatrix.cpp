#include "SimilarityMatrix.h"

/**
 * Constructor.
 */
SimilarityMatrix::SimilarityMatrix(EMatrix *ematrix, int quiet, char ** method,
    int num_methods, char * th_method, int x_coord, int y_coord,
    char * gene1, char * gene2, float th) {

  // Set some default values;
  this->ematrix = ematrix;
  this->quiet = quiet;
  this->method = method;
  this->num_methods = num_methods;
  this->th_method = th_method,
  this->x_coord = x_coord;
  this->y_coord = y_coord;
  this->gene1 = gene1;
  this->gene2 = gene2;
  this->th = th;

  // Find the index of the th_method in the methods array
  for (int i = 0; i < this->num_methods; i++) {
    if (strcmp(method[i], th_method) == 0) {
      this->th_method_index = i;
    }
  }

  if ((gene1 && !gene2) || (!gene1 && gene2)) {
    fprintf(stderr, "You must provide both gene1 and gene2 options.\n");
    exit(-1);
  }

  // if the user supplied gene
  if (gene1 && gene2) {
    // Make sure the coordinates are positive integers
    if (x_coord < 1) {
      fprintf(stderr, "Could not find gene %s in the genes list file\n", gene1);
      exit(-1);
    }
    if (y_coord < 1) {
      fprintf(stderr, "Could not find gene %s in the genes list file\n", gene2);
      exit(-1);
    }

    // If the user provided gene names then map those to the numeric IDs.
    char ** genes = ematrix->getGenes();
    int num_genes = ematrix->getNumGenes();
    int i = 0;
    for (i = 0; i < num_genes; i++) {
      if (strcmp(genes[i], gene1) == 0) {
        x_coord = i + 1;
      }
      if (strcmp(genes[i], gene2) == 0) {
         y_coord = i + 1;
      }
    }
  }

  // Make sure we have a positive integer for the x and y coordinates.
  if ((x_coord >= 1 &&  y_coord < 1) ||
      (x_coord < 1  &&  y_coord >= 1)) {
    fprintf(stderr, "Please provide a positive integer for both the x and y coordinates (-x and -y options)\n");
    exit(-1);
  }
}

/**
 * Destructor.
 */
SimilarityMatrix::~SimilarityMatrix() {

}
