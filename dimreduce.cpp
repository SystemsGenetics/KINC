#include "dimreduce.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

/**
 * DRArgs constructor.
 */
DRArgs::DRArgs(int argc, char *argv[]) {
  // Set defaults for all the private members.
  omit_na = 0;
  headers = 0;
  rows = 0;
  cols = 0;
  do_log10 = 0;
  do_log2 = 0;
  do_log = 0;
  min_obs = 30;
  msc_bw1 = 0.75;
  msc_bw2 = 0.9;
  infilename = NULL;
  na_val = NULL;

  // Initialize the 'func' parameter.
  strcpy(func, "none");

  // The value returned by getopt_long.
  int c;

  // Loop through the incoming arguments until the
  // getopt_long function returns -1. Then we break out of the loop.
  while(1) {
    int option_index = 0;

    // Specify the long options. The values returned are specified to be the
    // short options which are then handled by the case statement below.
    const struct option long_options[] = {
      {"omit_na",  no_argument,       &omit_na,  1 },
      {"headers",  no_argument,       &headers,  1 },
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
        infilename = optarg;
        break;
      case 'r':
        rows = atoi(optarg);
        break;
      case 'c':
        cols = atoi(optarg);
        break;
      case 'm':
        strcpy(method, optarg);
        break;
      case 'o':
        min_obs = atoi(optarg);
        break;
      case 'n':
        na_val = optarg;
        break;
      case 'f':
        strcpy(func, optarg);
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
  if (!infilename) {
    fprintf(stderr, "Please provide an expression matrix (--ematrix option).\n");
    exit(-1);
  }

  if (!method) {
    fprintf(stderr,"Please provide the method (--method option).\n");
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

  // make sure the input file exists
  if (access(infilename, F_OK) == -1) {
    fprintf(stderr, "Error: The input file does not exists or is not readable.\n");
    exit(-1);
  }

  // make sure the method is valid
  if (strcmp(method, "pc") != 0 &&
      strcmp(method, "mi") != 0 &&
      strcmp(method, "sc") != 0 ) {
    fprintf(stderr,"Error: The method (--method option) must either be 'pc', 'sc' or 'mi'.\n");
    exit(-1);
  }

  if (omit_na && !na_val) {
    fprintf(stderr, "Error: The missing value string should be provided (--na_val option).\n");
    exit(-1);
  }

  // TODO: make sure means shift bandwidth arguments are numeric between 1 and 0

  if (headers) {
    printf("  Skipping header lines\n");
  }
  printf("  Performing transformation: %s \n", func);
  printf("  Required observations: %d\n", min_obs);
  printf("  Using method: '%s'\n", method);
  if (omit_na) {
    printf("  Missing values are: '%s'\n", na_val);
  }

  // remove the path and extension from the filename
  char * filename = basename(infilename);
  strcpy(fileprefix, filename);
  char * p = rindex(fileprefix, '.');
  if (p) {
    p[0] = 0;
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

  // if we have a header line in the input file then
  // subtract one from the number of rows
  if (headers) {
    rows--;
  }

}

/**
 * DRArgs destructor.
 */
DRArgs::~DRArgs() {

}
/**
 *
 */
int do_dimreduce(int argc, char *argv[], int mpi_id, int mpi_num_procs) {

  // variables used for timing of the software
  time_t start_time = time(0);
  time_t now;

  DRArgs * params = new DRArgs(argc, argv);

  // Load the input expression matrix.
  EMatrix * ematrix = new EMatrix(params->getInfileName(), params->getNumRows(),
      params->getNumCols(), params->getHasHeaders(), params->getOmitNA(),
      params->getNAval());
  int i, j;

  // Open the clustering file for writing. Use the MPI ID in the filename.
//  FILE ** fps = open_output_files(params, mpi_id);

  // Calculate the total number of comparisons and how many will be
  // performed by this process. We subtract 1 from the first params->rows
  // because we do not calculate the diagonal.
  long long int num_rows = params->getNumRows() - 1;
  long long int total_comps  = num_rows * (num_rows + 1) / 2;
  long long int comps_per_process = total_comps / mpi_num_procs;
  long long int comp_start = mpi_id * comps_per_process;
  long long int comp_stop = mpi_id * comps_per_process + comps_per_process;

  // If this is the last process and there are some remainder comparisons
  // then we need to add them to the stop
  if (mpi_id + 1 == mpi_num_procs) {
    comp_stop = total_comps;
  }
  printf("%d. Performing %lld comparisons\n", mpi_id + 1, comp_stop - comp_start);
  fflush(stdout);

  // Perform the pair-wise clustering.
  int n_comps = 0;
  int my_comps = 0;
  for (i = 0; i < num_rows; i++) {
    /*if (i == 50) {
      break;
    }*/
    for (j = 0; j < num_rows; j++) {

      // We only need to calculate clusters in the lower triangle of the
      // full pair-wise matrix
      if (j >= i) {
        continue;
      }

      // If this computation is not meant for this process, then skip it.
      if (n_comps < comp_start || n_comps >= comp_stop) {
        n_comps++;
        continue;
      }

      if (n_comps % 100 == 0) {
        // Get the amount of memory used.
        statm_t * memory = memory_get_usage();

        // Calculate the number of days left.
        now = time(0);
        double seconds_passed = now - start_time;
        double minutes_passed = (seconds_passed) / 60.0;
        double percent_complete = (my_comps / (float) (comp_stop - comp_start)) * 100;
        double comps_per_minute = my_comps / minutes_passed;
        double total_time = (comp_stop - comp_start) / comps_per_minute;
        double minutes_left = total_time - minutes_passed;
        double hours_left = minutes_left / 60;
        double days_left = hours_left / 24;

        // Write progress report.
        printf("%d. Complete: %.4f%%. Mem: %ldb. Remaining: %.2fh; %.2fd. Coords: %d, %d.        \n",
          mpi_id + 1,
          (float) percent_complete,
          memory->size,
          (float) hours_left,
          (float) days_left,
          i,
          j
        );
        free(memory);
      }
      n_comps++;
      my_comps++;

      double *a = ematrix->getRow(j);
      double *b = ematrix->getRow(i);

      // Initialize the arrays that will be used for containing a and b
      // but with missing values removed.
      int cols = params->getNumCols();
      double * a2 = (double *) malloc(sizeof(double) * cols);
      double * b2 = (double *) malloc(sizeof(double) * cols);

      // The kept variable is an array of zeros and ones and is the size of the
      // number of samples.  A zero indicates the sample is not included in
      // the pair-wise comparison and a one indicates it is.
      int * kept = (int *) malloc(sizeof(int) * cols);
      int n2;

      // Remove any missing values before calculating Royston's H test.
      // This call will return a2 and b2 which are the same as a and b but
      // with missing values removed. n2 gets set to the size of a2 and b2.
      remove_missing_paired(a, b, cols, a2, b2, &n2, kept);

      /*// Perform the clustering if we have enough samples.
      if (n2 > 0) {

        // Perform the clustering.
        PairWiseClusters * clusters = clustering(a2, i, b2, j, n2, ematrix, params, 0.075, 0);

        // update the clusters to include zeros for any samples with a missing value.
        update_pairwise_cluster_samples(kept, params->cols, clusters);

        // Write the clusters to the file.
        write_pairwise_cluster_samples(clusters, fps);

        // clean up the memory.
        free_pairwise_cluster_list(clusters);
        free(kept);
      }
      // One of these genes has no observations so, just write out an
      // empty cluster.
      else {
        // Create an empty clusters set
        PairWiseClusters * newc = new_pairwise_cluster_list();
        newc->gene1 = i;
        newc->gene2 = j;
        newc->num_samples = params->cols;
        newc->samples = kept;
        newc->next = NULL;
        newc->cluster_size = 0;
        newc->pcc = NAN;
        write_pairwise_cluster_samples(newc, fps);
        free_pairwise_cluster_list(newc);
      }*/

      // Release the memory for a2 and b2.
      free(a2);
      free(b2);

    }
  }
  free(ematrix);
  free(params);
  return 1;
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
  printf("  --method|-m  The correlation method to use. Supported methods include\n");
  printf("                 Pearson's correlation ('pc'), Spearman's rank correlation ('sc')\n");
  printf("                 and Mutual Information ('mi'). Provide either 'pc', 'sc', or\n");
  printf("                 'mi' as values respectively.\n");
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
