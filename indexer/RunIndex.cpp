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
  printf("\n");
  printf("The list of required options:\n");
  printf("  --outdir|-o     The KINC output directory from a prevous run.\n");
  printf("  --samples|-g    The number of samples in the original ematrix.\n");
  printf("                  This corresponds to the --cols argument of the\n");
  printf("                  similarity program of KINC\n");
  printf("\n");
  printf("For Help:\n");
  printf("  --help|-h       Print these usage instructions\n");
  printf("\n");
}

/**
 * The function to call when running the 'index' command.
 */
RunIndex::RunIndex(int argc, char *argv[]) {
  nsamples = 0;

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
      {"outdir",       required_argument, 0,  'o' },
      // Expression matrix options.
      {"rows",         required_argument, 0,  'r' },
      {"cols",         required_argument, 0,  'c' },
      {"headers",      no_argument,       &headers,  1 },
      {"omit_na",      no_argument,       &omit_na,  1 },
      {"func",         required_argument, 0,  'f' },
      {"na_val",       required_argument, 0,  'n' },
      {"ematrix",      required_argument, 0,  'e' },

      // Last element required to be all zeros.
      {0, 0, 0,  0 }
     };

     // get the next option
     c = getopt_long(argc, argv, "o:h", long_options, &option_index);

     // if the index is -1 then we have reached the end of the options list
     // and we break out of the while loop
     if (c == -1) {
       break;
     }

     // handle the options
     switch (c) {
       case 0:
         break;
       case 'o':
         outdir = optarg;
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
         exit(-1);
     }
   }

  // Make sure an out file directory is provided
  if (!outdir) {
    fprintf(stderr, "Please provide the KINC output directory from a previous run (--outdir option).\n");
    exit(-1);
  }

  // Make sure the output directory exists.
  struct stat st = {0};
  if (stat(outdir, &st) == -1) {
    fprintf(stderr, "The specified output directory ,'%s', is missing. Please check the value of the --outdir argument.\n", outdir);
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

  // Retrieve the data from the EMatrix file.
  ematrix = new EMatrix(infilename, rows, cols, headers, omit_na, na_val, func);

  nsamples = ematrix->getNumSamples();

}
/**
 * Destructor
 */
RunIndex::~RunIndex() {
  //CLuceneIndexer indexer(outdir);
  SQLiteIndexer indexer(ematrix, outdir);
  indexer.run(nsamples);
}

/**
 * Performs the indexing.
 */
void RunIndex::execute() {

}

