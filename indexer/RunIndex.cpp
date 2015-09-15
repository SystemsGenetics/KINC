#include "../indexer/RunIndex.h"

/**
 * Prints the command-line usage instructions for the similarity command
 */
void RunIndex::printUsage() {
  printf("\n");
  printf("Usage: ./kinc index [options]\n");
  printf("The list of required options:\n");
  printf("  --ematrix|-e     The file name that contains the expression matrix.\n");
  printf("                   The rows must be genes or probe sets and columns are samples\n");
  printf("  --rows|-r        The number of lines in the ematrix file including the header\n");
  printf("                   row if it exists\n");
  printf("  --cols|-c        The number of columns in the input file\n");
  printf("  --method|-m      The correlation method to use. Supported methods include\n");
  printf("                   Pearson's correlation ('pc'), Spearman's rank correlation ('sc')\n");
  printf("                   and Mutual Information ('mi'). Provide either 'pc', 'sc', or\n");
  printf("                   'mi' as values respectively.\n");
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
  printf("Optional Arguments:\n");
  printf("  --job_index|-i   By default, indexing proceeds as a single job and processes\n");
  printf("                   one output directory at a time, however to index\n");
  printf("                   a single directory or to index in parallel, this \n");
  printf("                   option can be provided where the index is a numeric value \n");
  printf("                   between 0 and 101. If 101 is provided then the 'nan'\n");
  printf("                   directory is indexed.\n");
  printf("  --job_start|-s   This argument specifies the start job for indexing. Valid values\n");
  printf("                   are integers between 0 and 101.  The indexer will begin indexing\n");
  printf("                   the output directory with the given ID and proceed thereafter\n");
  printf("                   in descending order. This argument cannot be used at the same\n");
  printf("                   time as the job_index argument.\n");
  printf("\n");
  printf("For Help:\n");
  printf("  --help|-h        Print these usage instructions\n");
  printf("\n");
}

/**
 * The function to call when running the 'index' command.
 */
RunIndex::RunIndex(int argc, char *argv[]) {
  nsamples = 0;
  // Set a default value of -2 for the job_index. This means that
  // the indexing should proceed as a single job.
  job_index = -2;
  // Set the default job start to -2.  this means that no start is specified
  // by the user.
  job_start = -2;
  // The correlation method.
  method = NULL;

  // loop through the incoming arguments until the
  // getopt_long function returns -1. Then we break out of the loop
  // The value returned by getopt_long.
  int c;
  while(1) {
    int option_index = 0;

    // specify the long options. The values returned are specified to be the
    // short options which are then handled by the case statement below
    static struct option long_options[] = {
      {"help",         no_argument,       0,  'h' },
      {"method",       required_argument, 0,  'm' },
      // Expression matrix options.
      {"rows",         required_argument, 0,  'r' },
      {"cols",         required_argument, 0,  'c' },
      {"headers",      no_argument,       &headers,  1 },
      {"omit_na",      no_argument,       &omit_na,  1 },
      {"func",         required_argument, 0,  'f' },
      {"na_val",       required_argument, 0,  'n' },
      {"ematrix",      required_argument, 0,  'e' },
      // Index options.
      {"job_index",    required_argument, 0,  'i' },
      {"job_start",    required_argument, 0,  's' },
      // Last element required to be all zeros.
      {0, 0, 0,  0 }
     };

     // get the next option
     c = getopt_long(argc, argv, "o:m:r:c:f:n:e:i:s:h", long_options, &option_index);

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
         method = optarg;
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
       case 'i':
         job_index = atoi(optarg);
         break;
       case 's':
         job_start = atoi(optarg);
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
         exit(-1);
    }
  }
  if (!method) {
    fprintf(stderr, "Please provide the method (--method option) used to construct the similarity matrix.\n");
    exit(-1);
  }
  // make sure the method is valid
  if (strcmp(method, "pc") != 0 &&
      strcmp(method, "mi") != 0 &&
      strcmp(method, "sc") != 0 ) {
    fprintf(stderr,"Error: The method (--method option) must either be 'pc', 'sc' or 'mi'.\n");
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
  if (job_index < -2 || job_index > 101) {
    fprintf(stderr, "Please provide a job index that is between 0 and 101 (--job_index option).\n");
    exit(-1);
  }

  if (job_index > -2 && job_start > -2) {
    fprintf(stderr, "The job_index and job_start arguments cannot both be used at the same time).\n");
    exit(-1);
  }

  // Retrieve the data from the EMatrix file.
  ematrix = new EMatrix(infilename, rows, cols, headers, omit_na, na_val, func);

  nsamples = ematrix->getNumSamples();

  // Make sure the output directory exists.
  struct stat st = {0};
  sprintf(indexdir, "./clusters-%s", method);
  if (stat(indexdir, &st) == -1) {
    fprintf(stderr, "The specified indexes directory ,'%s', is missing. Please check the value of the --indexes argument.\n", indexdir);
    exit(-1);
  }

}
/**
 * Destructor
 */
RunIndex::~RunIndex() {
  //CLuceneIndexer indexer(indexdir);
  SQLiteIndexer indexer(ematrix, indexdir);
  indexer.run(nsamples, job_index, job_start);
}

/**
 * Performs the indexing.
 */
void RunIndex::execute() {

}

