#include "dimreduce.h"

int do_dimreduce(int argc, char *argv[]) {

  //time_t start_time, end_time;  // variables used for timing of the software
  int c;                        // the value returned by getopt_long

  static CCMParameters params;

  // initialize some of the program parameters
  params.perf = 1;
  params.omit_na = 0;
  params.headers = 0;
  params.rows = 0;
  params.cols = 0;
  params.do_log10 = 0;
  params.do_log2 = 0;
  params.do_log = 0;
  params.min_obs = 30;

  strcpy(params.func, "none");

  // loop through the incoming arguments until the
  // getopt_long function returns -1. Then we break out of the loop
  while(1) {
    int option_index = 0;

    // specify the long options. The values returned are specified to be the
    // short options which are then handled by the case statement below
    static struct option long_options[] = {
      {"omit_na",  no_argument,       &params.omit_na,  1 },
      {"perf",     no_argument,       &params.perf,     1 },
      {"headers",  no_argument,       &params.headers,  1 },
      {"help",     no_argument,       0,  'h' },
      {"ematrix",  required_argument, 0,  'e' },
      {"method",   required_argument, 0,  'm' },
      {"rows",     required_argument, 0,  'r' },
      {"cols",     required_argument, 0,  'c' },
      {"min_obs",  required_argument, 0,  'o' },
      {"func",     required_argument, 0,  'f' },
      {"na_val",   required_argument, 0,  'n' },
      {"mi_bins",  required_argument, 0,  'b' },
      {"mi_degree",required_argument, 0,  'd' },
      {0, 0, 0,  0 }  // last element required to be all zeros
    };
    // get the next option
    c = getopt_long(argc, argv, "e:r:c:m:o:n:f:h", long_options, &option_index);

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
        params.infilename = optarg;
        break;
      case 'r':
        params.rows = atoi(optarg);
        break;
      case 'c':
        params.cols = atoi(optarg);
        break;
      case 'o':
        params.min_obs = atoi(optarg);
        break;
      case 'n':
        params.na_val = optarg;
        break;
      case 'f':
        strcpy(params.func, optarg);
        break;
      case 'h':
        print_dimreduce_usage();
        exit(-1);
        break;
      case '?':
        exit(-1);
        break;
      case ':':
        print_dimreduce_usage();
        exit(-1);
        break;
      default:
        print_dimreduce_usage();
    }
  }

  // make sure the required arguments are set and appropriate
  if (!params.infilename) {
    fprintf(stderr,"Please provide an expression matrix (--ematrix option).\n");
    exit(-1);
  }
  // make sure we have a positive integer for the rows and columns of the matrix
  if (params.rows < 0 || params.rows == 0) {
    fprintf(stderr,"Please provide a positive integer value for the number of rows in the \n");
    fprintf(stderr,"expression matrix (--rows option).\n");
    exit(-1);
  }
  if (params.cols < 0 || params.cols == 0) {
    fprintf(stderr,"Please provide a positive integer value for the number of columns in\n");
    fprintf(stderr,"the expression matrix (--cols option).\n");
    exit(-1);
  }

  // make sure the input file exists
  if (access(params.infilename, F_OK) == -1) {
    fprintf(stderr,"Error: The input file does not exists or is not readable.\n");
    exit(-1);
  }

  if (params.omit_na && !params.na_val) {
    fprintf(stderr,"Error: The missing value string should be provided (--na_val option).\n");
    exit(-1);
  }

  if (params.headers) {
    printf("  Skipping header lines\n");
  }
  printf("  Performing transformation: %s \n", params.func);
  printf("  Required observations: %d\n", params.min_obs);
  if (params.omit_na) {
    printf("  Missing values are: '%s'\n", params.na_val);
  }

  // remove the path and extension from the filename
  char * filename = basename(params.infilename);
  strcpy(params.fileprefix, filename);
  char * p = rindex(params.fileprefix, '.');
  if (p) {
    p[0] = 0;
  }

  // if the function is log2 then set the do_log2 flag
  if (strcmp(params.func, "log10") == 0) {
    params.do_log10 = 1;
  }
  if (strcmp(params.func, "log2") == 0) {
    params.do_log2 = 1;
  }
  if (strcmp(params.func, "log") == 0) {
    params.do_log = 1;
  }

  // if we have a header line in the input file then
  // subtract one from the number of rows
  if (params.headers) {
    params.rows--;
  }

  EMatrix ematrix = load_ematrix(params);
  double ** data = ematrix.data;
  int total_comps;  // the total number of pair-wise comparisons to be made
  int n_comps;      // the number of comparisons completed during looping
  int i, j, k, l;

  total_comps = (params.rows * params.rows) / 2;
  for (i = 0; i < params.rows; i++) {
    for (j = 0; j < params.rows; j++) {

      n_comps++;
      if (n_comps % 1000 == 0) {
        printf("Percent complete: %.2f%%\r", (n_comps / (float) total_comps) * 100);
      }

      // We only need to calculate royson in the lower triangle of the
      // full pair-wise matrix
      if (j >= i) {
        continue;
      }
      double *a = data[i];
      double *b = data[j];

      // Initialize the arrays that will be used for containing a and b
      // but with missing values removed.
      double * a2 = (double *) malloc(sizeof(double) * params.cols);
      double * b2 = (double *) malloc(sizeof(double) * params.cols);
      int n2;

      // Remove any missing values before calculating Royston's H test.
      remove_missing_paired(a, b, params.cols, a2, b2, &n2);

      // Look for valid pair-wise comparison sets between these two genes.
      find_valid_comps(a2, i, b2, j, n2, ematrix, params);

      // Release the memory for a2 and b2.
      free(a2);
      free(b2);
    }
  }
  return 1;
}

/**
 *
 */
void find_valid_comps(double *a2, int x, double *b2, int y, int n2, EMatrix ematrix, CCMParameters params) {
  // Variables used for looping.
  int i, j, k, l;

  // Calculate Royston's H test for multivariate normality.
  double pv = royston2D(a2, b2, n2);
  printf("(%d, %d), pv: %e\n", i + 1, j + 1, pv);

  // If the Royston's H test has p-value < 0.05 which means
  // it does not appear to be multivariate normal, then do mean shift
  // clustering to cluster the measured points.
  if (pv < 0.05) {
    MeanShiftClusters clusters;

    // Calculate the clusters.
    clusters = meanshift2D(a2, b2, n2, 0.075);

    // Is there a cluster with the minimum observations?  If so,
    // then remove outliers.
    for(k = 0; k < clusters.num_clusters; k++) {
      if (clusters.sizes[k] >= params.min_obs) {

      }
    }

    // free the memory
    free(clusters.cluster_label);
    for(k = 0; k < clusters.num_clusters; k++) {
      for (l = 0; l < clusters.sizes[k]; l++) {
        free(clusters.clusters[k][l]);
      }
      free(clusters.clusters[k]);
    }
    free(clusters.sizes);
    free(clusters.clusters);
  }
  else {
    PairWiseSet pws;
    pws.gene1 = ematrix.samples[i];
    pws.gene2 = ematrix.samples[j];
    // we must know the sample names passed in for a2 and b2!!!!!
    write_clustering_file_line(pws)
  }
}

/**
 *
 */
void write_clustering_file_line(PairWiseSet pws) {

}

/**
 * Prints the command-line usage instructions for the similarity command
 */
void print_dimreduce_usage() {
  printf("\n");
  printf("Usage: ./kinc dimreduce [options]\n");
  printf("The list of required options:\n");
  printf("  --ematrix|-e The file name that contains the expression matrix.\n");
  printf("                 The rows must be genes or probe sets and columns are samples\n");
  printf("  --rows|-r    The number of lines in the input file including the header\n");
  printf("                 column if it exists\n");
  printf("  --cols|-c    The number of columns in the input file minus the first\n");
  printf("                 column that contains gene names\n");
  printf("\n");
  printf("Optional:\n");
  printf("  --omit_na     Provide this flag to ignore missing values. Defaults to 0.\n");
  printf("  --na_val|-n   A string representing the missing values in the input file\n");
  printf("                  (e.g. NA or 0.000)\n");
  printf("  --min_obs|-o  The minimum number of observations (after missing values\n");
  printf("                  removed) that must be present to perform correlation.\n");
  printf("                  Default is 30.\n");
  printf("  --func|-f     A transformation function to apply to elements of the ematrix.\n");
  printf("                  Values include: log, log2 or log10. Default is to not perform\n");
  printf("                  any transformation.\n");
  printf("  --perf        Provide this flag to enable performance monitoring.\n");
  printf("  --headers     Provide this flag if the first line of the matrix contains\n");
  printf("                  headers.\n");
  printf("  --help|-h     Print these usage instructions\n");
  printf("\n");
}
