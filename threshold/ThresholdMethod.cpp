#include "ThresholdMethod.h"

/**
 * Prints the command-line usage instructions for the similarity command
 */
void print_threshold_usage() {
  printf("\n");
  printf("Usage: ./kinc threshold [options]\n");
  printf("The list of required options:\n");
  printf("  --ematrix|-e  The file name that contains the expression matrix.\n");
  printf("                The rows must be genes or probe sets and columns are samples\n");
  printf("  --method|-m   The correlation method to use. Supported methods include\n");
  printf("                Pearson's correlation ('pc'), Spearman's rank correlation ('sc')\n");
  printf("                and Mutual Information ('mi'). Provide either 'pc', 'sc', or\n");
  printf("                'mi' as values respectively.\n");
  printf("  --rows|-r     The number of lines in the input file including the header\n");
  printf("                column if it exists\n");
  printf("  --cols|-c     The number of columns in the input file minus the first\n");
  printf("                column that contains gene names\n");
  printf("\n");
  printf("Optional:\n");
  printf("  --th|-t       A decimal indicating the start threshold. For Pearson's.\n");
  printf("                Correlation (--method pc), the default is 0.99. For Mutual\n");
  printf("                information (--method mi), the default is the maximum MI value\n");
  printf("                in the similarity matrix\n");
  printf("  --step|-s     The threshold step size, to subtract at each iteration of RMT.\n");
  printf("                The default is 0.001\n");
  printf("  --chi|-i      The Chi-square test value which when encountered RMT will stop.\n");
  printf("                The algorithm will only stop if it first encounters a Chi-square\n");
  printf("                value of 99.607 (df = 60, p-value = 0.001).\n");
  printf("  --missing|-g  The total number of allowed missing values.  Each gene\n");
  printf("                comparision can potentially have a missing value in any\n");
  printf("                sample.  When the maximum number of missing values exceeds\n");
  printf("                this value the clusters are ignored.\n");
  printf("                If not provided, no limit is set.\n");
  printf("  --size|-z     The minimum cluster size (number of samples per cluster).\n");
  printf("                Default is 30\n");
  printf("  --headers     Provide this flag if the first line of the matrix contains\n");
  printf("                headers.\n");
  printf("  --help|-h     Print these usage instructions\n");
  printf("\n");
}

/**
 * DRArgs constructor.
 */
ThresholdMethod::ThresholdMethod(int argc, char *argv[]) {
  // initialize some of the program parameters
  max_missing = INFINITY;
  min_cluster_size = 30;

  strcpy(method, "sc");

  // The value returned by getopt_long.
  int c;

  // loop through the incoming arguments until the
  // getopt_long function returns -1. Then we break out of the loop
  while(1) {
    int option_index = 0;

    // specify the long options. The values returned are specified to be the
    // short options which are then handled by the case statement below
    static struct option long_options[] = {
      {"help",    no_argument,       0,  'h' },
      {"ematrix", required_argument, 0,  'e' },
      {"method",  required_argument, 0,  'm' },
      {"missing", required_argument, 0,  'g' },
      {"size",    required_argument, 0,  'z' },
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
      case 'm':
        strcpy(method, optarg);
        break;
      case 'z':
        min_cluster_size = atoi(optarg);
        break;
      case 'g':
        max_missing = atoi(optarg);
        break;
      case 'h':
        print_threshold_usage();
        exit(-1);
        break;
      case '?':
        exit(-1);
        break;
      case ':':
        print_threshold_usage();
        exit(-1);
        break;
      default:
        print_threshold_usage();
    }
  }

  if (max_missing < -1) {
    fprintf(stderr, "Please provide a positive integer value for maximum missing values\n");
    fprintf(stderr, "the expression matrix (--max_missing option).\n");
    exit(-1);
  }
  if (min_cluster_size < 0 || min_cluster_size == 0) {
    fprintf(stderr, "Please provide a positive integer value for the minimum cluster size in\n");
    fprintf(stderr, "the expression matrix (--min_cluster_size option).\n");
    exit(-1);
  }

  if (!method) {
    fprintf(stderr,"Please provide the method used to construct the similarity matrix (--method option).\n");
    exit(-1);
  }
  if (strcmp(method, "pc") != 0 &&
      strcmp(method, "sc") != 0 &&
      strcmp(method, "mi") != 0) {
    fprintf(stderr,"The method (--method option) must either be 'pc', 'sc' or 'mi'.\n");
    exit(-1);
  }

  // Load the input expression matrix.
  ematrix = new EMatrix(argc, argv);

  input_dir = (char *) malloc(sizeof(char) * strlen(ematrix->getInfileName()));
  if (strcmp(method, "mi") == 0) {
    strcpy(input_dir, "MI");
  }
  else if (strcmp(method, "pc") == 0) {
    strcpy(input_dir, "Pearson");
  }
  else if (strcmp(method, "sc") == 0) {
    strcpy(input_dir, "Spearman");
  }
}

/**
 * DRArgs destructor.
 */
ThresholdMethod::~ThresholdMethod() {
  free(input_dir);
}
