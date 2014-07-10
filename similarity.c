#include "similarity.h"

/**
 * The function to call when running the 'similarity' command.
 */
int do_similarity(int argc, char *argv[]) {

  time_t start_time, end_time;  // variables used for timing of the software
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
  strcpy(params.method, "pc");

  // loop through the incoming arguments until the
  // getopt_long function returns -1. Then we break out of the loop
  while(1) {
    int option_index = 0;

    // specify the long options. The values returned are specified to be the
    // short options which are then handled by the case statement below
    static struct option long_options[] = {
      {"omit_na", no_argument,       &params.omit_na,  1 },
      {"perf",    no_argument,       &params.perf,     1 },
      {"headers", no_argument,       &params.headers,  1 },
      {"ematrix", required_argument, 0,  'e' },
      {"method",  required_argument, 0,  'm' },
      {"rows",    required_argument, 0,  'r' },
      {"cols",    required_argument, 0,  'c' },
      {"min_obs", required_argument, 0,  'o' },
      {"func",    required_argument, 0,  'f' },
      {"na_val",  required_argument, 0,  'n' },
      {0, 0, 0,  0 }  // last element required to be all zeros
    };

    // get the next option
    c = getopt_long(argc, argv, "e:r:c:m:o:n:f:", long_options, &option_index);

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
      case 'm':
        strcpy(params.method, optarg);
        printf("  Using method: '%s'\n", params.method);
        break;
      case 'o':
        params.min_obs = atoi(optarg);
        printf("  Requiring %d observations to perform pairwise calculations\n", params.min_obs);
        break;
      case 'n':
        params.na_val = optarg;
        printf("  Excluding missing values marked as '%s'\n", params.na_val);
        break;
      case 'f':
        strcpy(params.func, optarg);
        printf("  Performing %s transformation on expression matrix\n", params.func);
        break;
      case '?':
        exit(-1);
        break;
      case ':':
        print_similarity_usage();
        exit(-1);
        break;
      default:
        print_similarity_usage();
    }
  }

  if (params.omit_na) {
    printf("  Excluding missing values\n");
  }
  if (params.headers) {
    printf("  Skipping header\n");
  }

  // make sure the required arguments are set and appropriate
  if (!params.infilename) {
    printf("Please provide an expression matrix (--ematrix option). Use the -h option for help.\n");
    exit(-1);
  }

  if (!params.method) {
    printf("Please provide the method (--method option). Use the -h option for help.\n");
    exit(-1);
  }

  // make sure we have a positive integer for the rows and columns of the matrix
  if (params.rows < 0 || params.rows == 0) {
    printf("Please provide a positive integer value for the number of rows in the \n");
    printf("expression matrix (--rows option). Use the -h option for help.\n");
    exit(-1);
  }
  if (params.cols < 0 || params.cols == 0) {
    printf("Please provide a positive integer value for the number of columns in\n");
    printf("the expression matrix (--cols option). Use the -h option for help.\n");
    exit(-1);
  }

  // make sure the input file exists
  if (access(params.infilename, F_OK) == -1) {
    printf("Error: The input file does not exists or is not readable.\n");
    exit(-1);
  }

  if (strcmp(params.method, "pc") != 0 && strcmp(params.method, "mi") != 0) {
    printf("Error: The method (--method option) must either be 'pc' or 'mi'. Use the -h option for help.\n");
    exit(-1);
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

  // if we have a header line in the input file then subtract one from the number of rows
  if (params.headers) {
    params.rows--;
  }

  // if output of the histogram is enabled then initialize it
  int * histogram;
  histogram = init_histogram(params);

  // if performance monitoring is enabled the set the start time
  if (params.perf) {
    time(&start_time);
  }

  double ** data = load_ematrix(params);

  if (strcmp(params.method, "pc") == 0) {
    calculate_pearson(params, data, histogram);
  }
  else if(strcmp(params.method, "mi") == 0) {
    calculate_MI(params, data, histogram);
  }

  // if performance monitoring is enabled then write the timing data
  if(params.perf) {
    time(&end_time);
    FILE * timingfile = fopen("timingdata.txt", "a");
    fprintf(timingfile, "CCM Runtime with %d x %d %s input dataset: %.2lf min\n", params.rows, params.cols, params.infilename, difftime(end_time, start_time)/60.0);
  }

  // if output of the histogram is enabled then print it
  print_histogram(params, histogram);

  printf("Done.\n");
  return 1;
}

/**
 * Reads in the expression matrix
 *
 * @param CCMParameters params
 *   An instance of the CCM parameters struct
 *
 * @return
 *   A pointer to a two-dimensional array of doubles
 */
double ** load_ematrix(CCMParameters params) {

  FILE *infile; // pointers to the input
  char gene_name[255];   // holds the gene or probe set name from the first column
  double ** data;
  int i, j;  // integers for looping

  // allocate the data array for storing the input expression matrix
  data = (double**) malloc(sizeof(double *) * params.rows);
  for (i = 0; i < params.rows; i++) {
    data[i] = malloc(sizeof(double) * params.cols);
  }

  // iterate through the lines of the expression matrix
  printf("Reading input file '%s'...\n", params.infilename);
  infile = fopen(params.infilename, "r");
  for (i = 0; i < params.rows; i++) {

    // skip over the header if one is provided
    if (i == 0 && params.headers) {
      printf("Skipping headers...\n");
      char element[1024];
      for (j = 0; j < params.cols; j++) {
        fscanf(infile, "%s", element);
      }
    }

    // the first entry on every line is a label string - read that in before the numerical data
    fscanf(infile, "%s", gene_name);

    // iterate over the columns of each row
    for (j = 0; j < params.cols; j++) {
      char element[50]; // used for reading in each element of the expression matrix
      if (fscanf(infile, "%s", element) != EOF) {
        // if this is a missing value and omission of missing values is enabled then
        // rewrite this value as MISSING_VALUE
        if (params.omit_na && strcmp(element, params.na_val) == 0) {
          data[i][j] = NAN;
        }
        else {
          if (params.do_log10) {
            data[i][j] = log10(atof(element));
          }
          else if (params.do_log2) {
            data[i][j] = log2(atof(element));
          }
          else if (params.do_log) {
            data[i][j] = log(atof(element));
          }
          else {
            data[i][j] = atof(element);
          }
        }
      }
      else {
        printf("Error: EOF reached early. Exiting.\n");
        exit(-1);
      }
    }
  }
  return data;
}

/**
 * Creates an integer array, initialized with all values to zero for
 * storing a histogram. The number of bins in the histogram is
 * specified by the HIST_BINS global variable
 *
 * @param CCMParameters params
 *   An instance of the CCM parameters struct
 *
 * @return
 *   A pointer to an array of integers
 */
int * init_histogram(CCMParameters params) {
  int * histogram = (int *) malloc(sizeof(int) * HIST_BINS + 1);

  int m; // integer used or looping

  // create and initialize the histogram array for the distribution of coefficients
  for (m = 0; m < HIST_BINS + 1; m++) {
    histogram[m] = 0;
  }
  return histogram;
}

/**
 * Prints the histogram to a file.
 *
 * @param CCMParameters params
 *   An instance of the CCM parameters struct
 *
 * @param int * histogram
 *   An integer array created by the init_histogram function and
 *   populated by calculate_MI or calculate_person.
 *
 */
void print_histogram(CCMParameters params, int * histogram) {
  int m; // integer used or looping
  char outfilename[50];  // the output file name


  // output the correlation histogram
  sprintf(outfilename, "%s.mi.corrhist.txt", params.fileprefix);
  FILE * outfile = fopen(outfilename, "w");
  for (m = 0; m < HIST_BINS; m++) {
    fprintf(outfile, "%lf\t%d\n", 1.0 * m / (HIST_BINS), histogram[m]);
  }
  fclose(outfile);
}

/**
 * Prints the command-line usage instructions for the similarity command
 */
void print_similarity_usage() {
  printf("\n");
  printf("Usage: ./RMTGeneNet similarity [options]\n");
  printf("The list of required options:\n");
  printf("  --ematrix|-e The file name that contains the expression matrix.\n");
  printf("                 The rows must be genes or probe sets and columns are samples\n");
  printf("  --rows|-r    The number of lines in the input file including the header\n");
  printf("                 column if it exists\n");
  printf("  --cols|-c    The number of columns in the input file minus the first\n");
  printf("                 column that contains gene names\n");
  printf("  --method|-m  The correlation method to use. Supported methods include\n");
  printf("                 Pearson's correlation and Mutual Information. Provide\n");
  printf("                 either 'pc' or mi' as values respectively.\n");
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
  printf("\n");
  printf("Note: similarity values are set to NaN if there weren't enough observations\n");
  printf("to perform the calculation.\n");
}

