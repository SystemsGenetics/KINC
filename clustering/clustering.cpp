#include "clustering.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2


// Set to 1 if Ctrl-C is pressed. Will allow for graceful termination
int Close_Now = 0;

// Define the function to be called when ctrl-c (SIGINT) signal is sent to process
void signal_callback_handler(int signum) {
   printf("Caught signal %d. Terminating program...\n",signum);
   // Cleanup and close up stuff here
   Close_Now = 1;

}

/**
 * Prints the command-line usage instructions for the similarity command
 */
void print_clustering_usage() {
  printf("\n");
  printf("Usage: ./kinc dimreduce [options]\n");
  printf("The list of required options:\n");
  printf("  --ematrix|-e The file name that contains the expression matrix.\n");
  printf("               The rows must be genes or probe sets and columns are samples\n");
  printf("  --rows|-r    The number of lines in the input file including the header\n");
  printf("               column if it exists\n");
  printf("  --cols|-c    The number of columns in the input file minus the first\n");
  printf("               column that contains gene names\n");
  printf("  --method|-m  The correlation method to use. Supported methods include\n");
  printf("               Pearson's correlation ('pc'), Spearman's rank correlation ('sc')\n");
  printf("               and Mutual Information ('mi'). Provide either 'pc', 'sc', or\n");
  printf("               'mi' as values respectively.\n");
  printf("\n");
  printf("Mixture Models Options:\n");
  printf("  --criterion|-t     The criterion model to use: BIC, ICL, NEC, CV, DCV.\n");
  printf("  --max_clusters|-l  The maximum number of clusters to look for.\n");
  printf("\n");
  printf("Optional:\n");
  printf("  --num_jobs   Dimensionality reduction using clustering is highly parallel,\n");
  printf("               thus the task can be split into multiple separate jobs. KINC \n");
  printf("               can be run multiple times each working on a different portion \n");
  printf("               of the total comparisions.  Use this argument to specify the \n");
  printf("               total number of separate jobs that will be executed in parallel.\n");
  printf("               KINC will not launch these. They must be launched manually, each\n");
  printf("               with a different job_index. Defaults to 1\n");
  printf("  --job_index  When the argument num_jobs is provided, this argument indicates which\n");
  printf("               job is being run. For example, if num_jobs is set to 10, then this\n");
  printf("               argument should be any number between 1 and 10. Defaults to 1.\n");
  printf("  --omit_na    Provide this flag to ignore missing values. Defaults to 0.\n");
  printf("  --na_val|-n  A string representing the missing values in the input file\n");
  printf("               (e.g. NA or 0.000)\n");
  printf("  --min_obs|-o The minimum number of observations (after missing values\n");
  printf("               removed) that must be present to perform correlation.\n");
  printf("               Default is 30.\n");
  printf("  --func|-f    A transformation function to apply to elements of the ematrix.\n");
  printf("               Values include: log, log2 or log10. Default is to not perform\n");
  printf("               any transformation.\n");
  printf("  --perf       Provide this flag to enable performance monitoring.\n");
  printf("  --headers    Provide this flag if the first line of the matrix contains\n");
  printf("               headers.\n");
  printf("  --help|-h    Print these usage instructions\n");
  printf("\n");
}

/**
 *
 */
//int do_dimreduce(int argc, char *argv[], int mpi_id, int mpi_num_procs) {
void do_pairwise_Clustering(int argc, char *argv[], EMatrix *ematrix) {

  MixtureModelPWClustering mixmodc = new MixtureModelPWClustering(argc, argv, ematrix);

  mixmodc->run();
}
/**
 * DRArgs constructor.
 */
PairWiseClustering::PairWiseClustering(int argc, char *argv[], EMatrix *ematrix) {
  // Set defaults for all the private members.
  min_obs = 30;
  num_jobs = 1;
  job_index = 1;
  this->ematrix = ematrix;

  // The value returned by getopt_long.
  int c;

  // Loop through the incoming arguments until the
  // getopt_long function returns -1. Then we break out of the loop.
  while(1) {
    int option_index = 0;

    // Specify the long options. The values returned are specified to be the
    // short options which are then handled by the case statement below.
    const struct option long_options[] = {
      {"help",         no_argument,       0,  'h' },
      {"method",       required_argument, 0,  'm' },
      {"num_jobs",     required_argument, 0,  'j' },
      {"job_index",    required_argument, 0,  'i' },
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
      case 'm':
        strcpy(method, optarg);
        break;
      case 'o':
        min_obs = atoi(optarg);
        break;
      case 'j':
        num_jobs = atoi(optarg);
        break;
      case 'i':
        job_index = atoi(optarg);
        break;
      case 'h':
        print_clustering_usage();
        exit(-1);
        break;
      case '?':
        exit(-1);
        break;
      case ':':
        print_clustering_usage();
        exit(-1);
        break;
      default:
        print_clustering_usage();
    }
  }

  if (!method) {
    fprintf(stderr,"Please provide the method (--method option).\n");
    exit(-1);
  }

  // make sure the method is valid
  if (strcmp(method, "pc") != 0 &&
      strcmp(method, "mi") != 0 &&
      strcmp(method, "sc") != 0 ) {
    fprintf(stderr,"Error: The method (--method option) must either be 'pc', 'sc' or 'mi'.\n");
    exit(-1);
  }

  // TODO: make sure means shift bandwidth arguments are numeric between 1 and 0

  // Load the input expression matrix.
  EMatrix * ematrix = new EMatrix(argc, argv);

  printf("  Required observations: %d\n", min_obs);
  printf("  Using method: '%s'\n", method);
}

/**
 * DRArgs destructor.
 */
PairWiseClustering::~PairWiseClustering() {
  delete ematrix;
}
