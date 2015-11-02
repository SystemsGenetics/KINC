#include "RunThreshold.h"

/**
 * Prints the command-line usage instructions for the similarity command
 */
void RunThreshold::printUsage() {
  printf("\n");
  printf("Usage: ./kinc threshold [options]\n");
  printf("The list of required options:\n");
  printf("  --ematrix|-e     The file name that contains the expression matrix.\n");
  printf("                   The rows must be genes or probe sets and columns are samples\n");
  printf("  --method|-m      The correlation methods used. Supported methods include\n");
  printf("                   Pearson's correlation ('pc'), Spearman's rank ('sc')\n");
  printf("                   and Mutual Information ('mi'). Provide the methods that\n");
  printf("                   were used in the 'similarity' step.  The methods must appear\n");
  printf("                   in the same order.");
  printf("  --th_method|-t   The method to be used for thresholding.  Only one method can\n");
  printf("                   be used for thresholding even if multiple methods are provided\n");
  printf("                   for the --method argument.\n");
  printf("  --rows|-r        The number of lines in the input file including the header\n");
  printf("                   column if it exists\n");
  printf("  --cols|-c        The number of columns in the input file minus the first\n");
  printf("                   column that contains gene names\n");
  printf("\n");
  printf("Optional expression matrix arguments:\n");
  printf("  --omit_na        Provide this flag to ignore missing values.\n");
  printf("  --na_val|-n      A string representing the missing values in the input file\n");
  printf("                   (e.g. NA or 0.000)\n");
  printf("  --func|-f        A transformation function to apply to elements of the ematrix.\n");
  printf("                   Values include: log, log2 or log10. Default is to not perform\n");
  printf("                   any transformation.\n");
  printf("  --headers        Provide this flag if the first line of the matrix contains\n");
  printf("                   headers.\n");
  printf("\n");
  printf("Optional RMT arguments:\n");
  printf("  --th|-t          A decimal indicating the start threshold. For Pearson's.\n");
  printf("                   and Spearman's, the default is 0.99. For Mutual\n");
  printf("                   information (--method mi), the default is the maximum MI value\n");
  printf("                   in the similarity matrix\n");
  printf("  --step|-s        The threshold step size, to subtract at each iteration of RMT.\n");
  printf("                   The default is 0.001\n");
  printf("  --chi|-i         The Chi-square test value which when encountered RMT will stop.\n");
  printf("                   The algorithm will only stop if it first encounters a Chi-square\n");
  printf("                   value of 99.607 (df = 60, p-value = 0.001).\n");
  printf("\n");
  printf("Optional clustered data arguments:\n");
  printf("  --clustering|-l  The type of clustering that was performed during construction\n");
  printf("                   of the similarity matrix (e.g. 'mixmod').\n");
  printf("  --missing|-g     The total number of allowed missing values.  Each gene\n");
  printf("                   comparision can potentially have a missing value in any\n");
  printf("                   sample.  When the maximum number of missing values exceeds\n");
  printf("                   this value the clusters are ignored.\n");
  printf("                   If not provided, no limit is set.\n");
  printf("  --min_csize|-z   The minimum cluster size (number of samples per cluster).\n");
  printf("                   Default is 30\n");
  printf("  --max_modes|-d   The maximum number of modes. If a pair-wise comparision\n");
  printf("                   has multiple modes (i.e. multiple clusters) then only clusters from those\n");
  printf("                   comparisions with modes equal to or less than the value specified are\n");
  printf("                   included. Default is 1.\n");
  printf("\n");
  printf("For Help:\n");
  printf("  --help|-h     Print these usage instructions\n");
  printf("\n");
}

/**
 *
 */
RunThreshold::RunThreshold(int argc, char *argv[]) {
  // initialize some of the program parameters
  max_missing = INFINITY;
  max_modes = 1;
  min_cluster_size = 30;

  thresholdStart = 0.99;
  thresholdStep  = 0.001;
  chiSoughtValue = 200;

  // Initialize the array of method names. We set it to 10 as max. We'll
  // most likely never have this many of similarity methods available.
  method = (char **) malloc(sizeof(char *) * 10);

  // The concatenated method.
  char * cmethod;

  // The value returned by getopt_long.
  int c;

  // loop through the incoming arguments until the
  // getopt_long function returns -1. Then we break out of the loop
  while(1) {
    int option_index = 0;

    // specify the long options. The values returned are specified to be the
    // short options which are then handled by the case statement below
    static struct option long_options[] = {
      {"help",         no_argument,       0,  'h' },
      {"method",       required_argument, 0,  'm' },
      {"th_method",    required_argument, 0,  'p' },
      {"missing",      required_argument, 0,  'g' },
      {"size",         required_argument, 0,  'z' },
      // Expression matrix options.
      {"rows",         required_argument, 0,  'r' },
      {"cols",         required_argument, 0,  'c' },
      {"headers",      no_argument,       &headers,  1 },
      {"omit_na",      no_argument,       &omit_na,  1 },
      {"func",         required_argument, 0,  'f' },
      {"na_val",       required_argument, 0,  'n' },
      {"ematrix",      required_argument, 0,  'e' },
      // Clustered options
      {"clustering",   required_argument, 0,  'l' },
      {"max_missing",  required_argument, 0,  'g' },
      {"min_csize",    required_argument, 0,  'z' },
      {"max_modes",    required_argument, 0,  'd' },
      // RMT Threshold options.
      {"chi",          required_argument, 0,  'i' },
      {"th",           required_argument, 0,  't' },
      {"step",         required_argument, 0,  's' },

      // Last element required to be all zeros.
      {0, 0, 0,  0 }
    };

    // get the next option
    c = getopt_long(argc, argv, "m:g:z:r:c:f:n:e:t:d:l:h", long_options, &option_index);

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
        cmethod = optarg;
        break;
      case 'p':
        th_method = optarg;
        break;
      // Clustering options.
      case 'l':
        clustering = optarg;
        break;
      case 'z':
        min_cluster_size = atoi(optarg);
        break;
      case 'g':
        max_missing = atoi(optarg);
        break;
      case 'd':
        max_modes = atoi(optarg);
        break;
      // RMT threshold options.
      case 't':
        thresholdStart = atof(optarg);
        break;
      case 'i':
        chiSoughtValue = atof(optarg);
        break;
      case 's':
        thresholdStep = atof(optarg);
        break;
      // Expression matrix options.
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
      // Help and catch-all options.
      case 'h':
        printUsage();
        exit(-1);
        break;
      case '?':
        exit(-1);
        break;
      case ':':
        printUsage();
        exit(-1);
        break;
      default:
        printUsage();
    }
  }

  // Make sure the similarity method is valid.
  if (!cmethod) {
    fprintf(stderr,"Please provide the method (--method option).\n");
    exit(-1);
  }
  parseMethods(cmethod);

  if (clustering) {
    // Make sure the clustering method is valid.
    if (strcmp(clustering, "mixmod") != 0) {
      fprintf(stderr,"Error: The clustering type (--clustering option) must be 'mixmod'. It is currently the only clustering method supported.\n");
          exit(-1);
    }

    if (max_missing < -1) {
      fprintf(stderr, "Please provide a positive integer value for maximum missing values (--max_missing option).\n");
      exit(-1);
    }
    if (max_modes < 1) {
      fprintf(stderr, "Please provide a positive integer greater than 1 for the maximum modes values (--max_modes option).\n");
      exit(-1);
    }
    if (min_cluster_size < 0 || min_cluster_size == 0) {
      fprintf(stderr, "Please provide a positive integer value for the minimum cluster size (--min_cluster_size option).\n");
      exit(-1);
    }
  }

  // Make sure we have a positive integer for the rows and columns of the matrix.
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

  // make sure the required arguments are set and appropriate
  if (!infilename) {
    fprintf(stderr,"Please provide an expression matrix (--ematrix option).\n");
    exit(-1);
  }

  // TODO: make sure the th_method is in the method array.

  // Load the input expression matrix.
  ematrix = new EMatrix(infilename, rows, cols, headers, omit_na, na_val, func);

  if (headers) {
    printf("  Skipping header lines\n");
  }
  printf("  Performing transformation: %s \n", func);
  if (omit_na) {
    printf("  Missing values are: '%s'\n", na_val);
  }
  for(int i = 0; i < this->num_methods; i++) {
    printf("  Expecting similarity methods: '%s'\n", method[i]);
  }
  printf("  Method for thresholding: %s\n", th_method);
  printf("  Start threshold: %f\n", thresholdStart);
  printf("  Stopping Chi-square %f\n", chiSoughtValue);
  printf("  Step per iteration: %f\n", thresholdStep);
}
/**
 *
 */
void RunThreshold::parseMethods(char * methods_str) {
  // Split the method into as many parts
  char * tmp;
  int i = 0;
  tmp = strstr(methods_str, ",");
  while (tmp) {
    // Get the method and make sure it's valid.
    method[i] = (char *) malloc(sizeof(char) * (tmp - methods_str + 1));
    strncpy(method[i], methods_str, (tmp - methods_str));
    methods_str = tmp + 1;
    tmp = strstr(methods_str, ",");
    i++;
  }
  // Get the last element of the methods_str.
  method[i] = (char *) malloc(sizeof(char) * strlen(methods_str) + 1);
  strcpy(method[i], methods_str);
  i++;

  this->num_methods = i;

  for (i = 0; i < this->num_methods; i++) {
    if (strcmp(method[i], "pc") != 0 &&
        strcmp(method[i], "mi") != 0 &&
        strcmp(method[i], "sc") != 0 ) {
      fprintf(stderr,"Error: The method (--method option) must contain only 'pc', 'sc' or 'mi'.\n");
      exit(-1);
    }
    // Make sure the method isn't specified more than once.
    for (int j = 0; j < i; j++) {
      if (strcmp(method[i], method[j]) == 0) {
        fprintf(stderr,"Error: You may only specify a similarity method once (--method option).\n");
        exit(-1);
      }
    }
  }
}
/**
 *
 */
RunThreshold::~RunThreshold() {
  delete ematrix;
}

/**
 *
 */
void RunThreshold::execute() {

  // Find the RMT threshold.
  RMTThreshold * rmt = new RMTThreshold(ematrix, method, num_methods, th_method, thresholdStart,
      thresholdStep, chiSoughtValue, clustering, min_cluster_size, max_missing,
      max_modes);
  rmt->findThreshold();
  printf("Done.\n");
}
