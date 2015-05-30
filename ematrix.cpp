#include "ematrix.h"

/**
 * Reads in the expression matrix
 *
 *
 * @return
 *   A pointer to a two-dimensional array of doubles
 */
EMatrix::EMatrix(int argc, char *argv[]) {

  // initialize some of the program parameters
  rows = 0;
  cols = 0;
  headers = 0;
  omit_na = 0;
  na_val = NULL;
  do_log10 = 0;
  do_log2 = 0;
  do_log = 0;

  // Initialize the 'func' parameter.
  strcpy(func, "none");

  // The value returned by getopt_long.
  int c;

  // loop through the incoming arguments until the
  // getopt_long function returns -1. Then we break out of the loop
  while(1) {
    int option_index = 0;

    // specify the long options. The values returned are specified to be the
    // short options which are then handled by the case statement below
    static struct option long_options[] = {
      {"rows",    required_argument, 0,  'r' },
      {"cols",    required_argument, 0,  'c' },
      {"help",    no_argument,       0,  'h' },
      {"headers", no_argument,       &headers,  1 },
      {"omit_na", no_argument,       &omit_na,  1 },
      {"func",    required_argument, 0,  'f' },
      {"na_val",  required_argument, 0,  'n' },
      {"ematrix", required_argument, 0,  'e' },
      {0, 0, 0,  0 }  // last element required to be all zeros
    };

    // get the next option
    c = getopt_long(argc, argv, "e:m:t:c:s:h", long_options, &option_index);

    // if the index is -1 then we have reached the end of the options list
    // and we break out of the while loop
    if (c == -1) {
      break;
    }

    // handle the options
    switch (c) {
      case 0:
        break;
      case 'e':
        infilename = optarg;
        break;
      case 'r':
        rows = atoi(optarg);
        break;
      case 'c':
        cols = atoi(optarg);
        break;
      case 'n':
        na_val = optarg;
        break;
      case 'f':
        strcpy(func, optarg);
        break;
    }
  }

  // make sure the required arguments are set and appropriate
  if (!infilename) {
    fprintf(stderr,"Please provide an expression matrix (--ematrix option).\n");
    exit(-1);
  }
  // make sure we have a positive integer for the rows and columns of the matrix
  if (rows < 0 || rows == 0) {
    fprintf(stderr, "Please provide a positive integer value for the number of rows in the \n");
    fprintf(stderr, "expression matrix (--rows option).\n");
    exit(-1);
  }
  if (cols < 0 || cols == 0) {
    fprintf(stderr, "Please provide a positive integer value for the number of columns in\n");
    fprintf(stderr, "the expression matrix (--cols option).\n");
    exit(-1);
  }

  if (omit_na && !na_val) {
    fprintf(stderr, "Error: The missing value string should be provided (--na_val option).\n");
    exit(-1);
  }
  // make sure the input file exists
  if (access(infilename, F_OK) == -1) {
    fprintf(stderr,"The input file does not exists or is not readable.\n");
    exit(-1);
  }

  // if the function is log2 then set the do_log2 flag
  if (strcmp(func, "log10") == 0) {
    do_log10 = 1;
  }
  if (strcmp(func, "log2") == 0) {
    do_log2 = 1;
  }
  if (strcmp(func, "log") == 0) {
    do_log = 1;
  }

  // Remove the path and extension from the filename.
  file_prefix = (char *) malloc(sizeof(char) * strlen(infilename));
  char * temp = basename((char *) infilename);
  strcpy(file_prefix, temp);
  char * p = rindex(file_prefix, '.');
  if (p) {
    p[0] = 0;
  }

  // Pointer to the input file.
  FILE *infile;
  // Integers for looping.
  int i, j;

  // TODO: make sure means shift bandwidth arguments are numeric between 1 and 0
  if (headers) {
    printf("  Skipping header lines\n");
  }

  printf("  Performing transformation: %s \n", func);
  if (omit_na) {
    printf("  Missing values are: '%s'\n", na_val);
  }

  // Initialize the genes and samples arrays.
  num_genes = rows;
  num_samples = cols;
  if (headers) {
    num_genes--;
    samples = (char **) malloc(sizeof(char *) * num_samples);
  }
  genes = (char **) malloc(sizeof(char *) * num_genes);

  // Allocate the data array for storing the input expression matrix.
  data = (double**) malloc(sizeof(double *) * num_genes);
  for (i = 0; i < num_genes; i++) {
    data[i] = (double *) malloc(sizeof(double) * num_samples);
  }

  // Iterate through the lines of the expression matrix
  printf("Reading input file '%s' with %d rows and %d columns...\n", infilename, rows, cols);
  infile = fopen(infilename, "r");
  int k = 0;
  for (i = 0; i < rows; i++) {

    // If a header is provided then get the sample names.
    if (i == 0 && headers) {
      printf("Reading sample names...\n");

      for (j = 0; j < num_samples; j++) {
        samples[j] = (char *) malloc(sizeof(char) * 255);
        fscanf(infile, "%s", samples[j]);
      }
      continue;
    }

    // The first entry on every line is a label string - read that in before the numerical data
    genes[k] = (char *) malloc(sizeof(char) * 255);
    fscanf(infile, "%s", genes[k]);

    // iterate over the columns of each row
    for (j = 0; j < num_samples; j++) {
      char element[50]; // used for reading in each element of the expression matrix
      if (fscanf(infile, "%s", element) != EOF) {
        // if this is a missing value and omission of missing values is enabled then
        // rewrite this value as MISSING_VALUE
        if (omit_na && strcmp(element, na_val) == 0) {
          data[k][j] = NAN;
        }
        else {
          // make sure the element is numeric
          if (!is_numeric(element)) {
            fprintf(stderr, "Error: value is not numeric: %s\n", element);
            exit(-1);
          }
          data[k][j] = atof(element);
        }
      }
      else {
        fprintf(stderr, "Error: EOF reached early. Exiting.\n");
        exit(-1);
      }
    }
    k++;
  }

  // Perform any transformations to the data requested by the user.
  if (do_log) {
    logTransform();
  }
  if (do_log2) {
    log2Transform();
  }
  if (do_log10) {
    log10Transform();
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

