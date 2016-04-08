#include "RunQuery.h"

/**
 * Prints the command-line usage instructions for the similarity command
 */
void RunQuery::printUsage() {
  printf("\n");
  printf("Usage: ./kinc query [options]\n");
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
  printf("Optional expression matrix arguments:\n");
  printf("  --omit_na        Provide this flag to ignore missing values.\n");
  printf("  --na_val|-n      A string representing the missing values in the input file\n");
  printf("                   (e.g. NA or 0.000)\n");
  printf("  --func|-f        A transformation function to apply to elements of the ematrix.\n");
  printf("                   Values include: log, log2 or log10. Default is to not perform\n");
  printf("                   any transformation.\n");
  printf("  --headers        Provide this flag if the first line of the matrix contains\n");
  printf("                   headers.\n");
  printf("Optional query arguments:\n");
  printf("  -x               Retrieve results for gene x.\n");
  printf("  -y               Retrieve results for gene y.\n");
  printf("  --gene1|-1       The name of the first gene in a singe\n");
  printf("                   pair-wise comparision.  Must be used with --gene2 option.\n");
  printf("  --gene2|-2       The name of the second gene in a singe\n");
  printf("                   pair-wise comparision.  Must be used with --gene1 option.\n");
  printf("  --score|-s       Retrieve results with a similarity score >= the\n");
  printf("                   absolute value of this argument.\n");
  printf("\n");
  printf("For Help:\n");
  printf("  --help|-h       Print these usage instructions\n");
  printf("\n");
}

/**
 * The function to call when running the 'index' command.
 */
RunQuery::RunQuery(int argc, char *argv[]) {
  headers = 0;
  rows = 0;
  cols = 0;
  x_coord = -1;
  y_coord = -1;
  gene1 = NULL;
  gene2 = NULL;
  score = 0;
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
      // Index query options.
      {"indexdir",     required_argument, 0,  'i' },
      {"score",        required_argument, 0,  's' },
      {"x",            required_argument, 0,  'x' },
      {"y",            required_argument, 0,  'y' },
      {"gene1",        required_argument, 0,  '1' },
      {"gene2",        required_argument, 0,  '2' },
      // Last element required to be all zeros.
      {0, 0, 0,  0 }
     };
    delete ematrix;

     // get the next option
     c = getopt_long(argc, argv, "m:r:c:f:n:e:i:x:y:s:1:2:h", long_options, &option_index);

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
       case 'x':
         x_coord = atoi(optarg);
         break;
       case 'y':
         y_coord = atoi(optarg);
         break;
       case '1':
         gene1 = optarg;
         break;
       case '2':
         gene2 = optarg;
         break;
       case 's':
         score = atof(optarg);
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

  // Load the input expression matrix.
  ematrix = new EMatrix(infilename, rows, cols, headers, omit_na, na_val, func);

  if ((gene1 && !gene2) || (!gene1 && gene2)) {
    fprintf(stderr, "You must provide both gene1 and gene2 options.\n");
    exit(-1);
  }
  // if the user supplied gene
  if (gene1 && gene2) {
    x_coord = ematrix->getGeneCoord(gene1);
    y_coord = ematrix->getGeneCoord(gene2);

    // Make sure the coordinates are positive integers
    if (x_coord < 0) {
      fprintf(stderr, "Could not find gene %s in the genes list file\n", gene1);
      exit(-1);
    }
    if (y_coord < 0) {
      fprintf(stderr, "Could not find gene %s in the genes list file\n", gene2);
      exit(-1);
    }
  }

  // Make sure we have a positive integer for the x and y coordinates.
  if ((x_coord >= 1 &&  y_coord < 1) ||
      (x_coord < 1  &&  y_coord >= 1)) {
    fprintf(stderr, "Please provide a positive integer greater than 1 for both the x and y coordinates (-x and -y options)\n");
    exit(-1);
  }

  if (score < 0 || score > 1) {
    fprintf(stderr, "Please provide a score between 0 and 1 (--s option).\n");
    exit(-1);
  }

  // Make sure the output directory exists.
//  struct stat st = {0};
//  sprintf(indexdir, "./clusters-%s", method);
//  if (stat(indexdir, &st) == -1) {
//    fprintf(stderr, "The specified indexes directory ,'%s', is missing. Please check the value of the --indexes argument.\n", indexdir);
//    exit(-1);
//  }
}
/**
 * Destructor
 */
RunQuery::~RunQuery() {

}

/**
 * Performs the indexing.
 */
void RunQuery::execute() {
//  CLuceneQuery query(indexdir);
//  SQLiteQuery query(indexdir, ematrix);
//  query.run(x_coord, y_coord, score);
}
