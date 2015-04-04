#include "ematrix.h"

/**
 * Reads in the expression matrix
 *
 * @param CCMParameters params
 *   An instance of the CCM parameters struct
 *
 * @return
 *   A pointer to a two-dimensional array of doubles
 */
EMatrix::EMatrix(char * infilename, int rows, int cols, int headers,
    int omit_na, char * na_val) {

  // Pointer to the input file.
  FILE *infile;
  // Integers for looping.
  int i, j;

  // Initialize the EMatrix struct
  num_genes = rows;
  num_samples = cols;
  if (headers) {
    num_genes--;
    samples = (char **) malloc(sizeof(char *) * num_samples);
  }
  genes = (char **) malloc(sizeof(char *) * num_genes);

  // Allocate the data array for storing the input expression matrix.
  data = (double**) malloc(sizeof(double *) * rows);
  for (i = 0; i < rows; i++) {
    data[i] = (double *) malloc(sizeof(double) * cols);
  }

  // iterate through the lines of the expression matrix
  printf("Reading input file '%s'...\n", infilename);
  infile = fopen(infilename, "r");
  for (i = 0; i < rows; i++) {

    // skip over the header if one is provided
    if (i == 0 && headers) {
      printf("Skipping headers...\n");

      for (j = 0; j < cols; j++) {
        samples[j] = (char *) malloc(sizeof(char) * 255);
        fscanf(infile, "%s", samples[j]);
      }
    }

    // the first entry on every line is a label string - read that in before the numerical data
    genes[i] = (char *) malloc(sizeof(char) * 255);
    fscanf(infile, "%s", genes[i]);

    // iterate over the columns of each row
    for (j = 0; j < cols; j++) {
      char element[50]; // used for reading in each element of the expression matrix
      if (fscanf(infile, "%s", element) != EOF) {
        // if this is a missing value and omission of missing values is enabled then
        // rewrite this value as MISSING_VALUE
        if (omit_na && strcmp(element, na_val) == 0) {
          data[i][j] = NAN;
        }
        else {
          // make sure the element is numeric
          if (!is_numeric(element)) {
            fprintf(stderr, "Error: value is not numeric: %s\n", element);
            exit(-1);
          }
          data[i][j] = atof(element);
        }
      }
      else {
        fprintf(stderr, "Error: EOF reached early. Exiting.\n");
        exit(-1);
      }
    }
  }
}
/**
 * Frees the memory associated with an EMatrix object.
 *
 * @param EMatrix ematrix
 *   An instance of the EMatrix struct.
 */
EMatrix::~EMatrix() {
  int i;
  for (i = 0; i < num_samples; i++) {
    free(samples[i]);
  }
  for (i = 0; i < num_genes; i++) {
    free(data[i]);
    free(genes[i]);
  }
  free(data);
  free(genes);
  free(samples);
}

/**
 *
 */
void EMatrix::logTransform() {
  for (int i = 0; i < num_genes; i++) {
    for (int j = 0; j < num_samples; j++) {
      data[i][j] = log(data[i][j]);
    }
  }
}
/**
 *
 */
void EMatrix::log2Transform() {
  for (int i = 0; i < num_genes; i++) {
    for (int j = 0; j < num_samples; j++) {
      data[i][j] = log2(data[i][j]);
    }
  }
}
/**
 *
 */
void EMatrix::log10Transform(){
  for (int i = 0; i < num_genes; i++) {
    for (int j = 0; j < num_samples; j++) {
      data[i][j] = log10(data[i][j]);
    }
  }
}

