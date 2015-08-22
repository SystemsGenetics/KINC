#include "RunQuery.h"

/**
 * Prints the command-line usage instructions for the similarity command
 */
void RunQuery::printUsage() {
  printf("\n");
  printf("Usage: ./kinc query [options]\n");
  printf("The list of required options:\n");
  printf("  --indexdir|-i    The KINC directory where the indexed files are housed.\n");
  printf("  --outfile|-o     The file where results should be written.\n");
  printf("  -x               Retrieve results for gene x.\n");
  printf("  -y               Retrieve results for gene y.\n");
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
  x_coord = 0;
  y_coord = 0;
  score = 0.92;

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
      {"indexdir",     required_argument, 0,  'i' },
      {"outfile",      required_argument, 0,  'o' },
      {"score",        required_argument, 0,  's' },
      {"x",            required_argument, 0,  'x' },
      {"y",            required_argument, 0,  'y' },
      // Last element required to be all zeros.
      {0, 0, 0,  0 }
     };

     // get the next option
     c = getopt_long(argc, argv, "i:x:y:s:h", long_options, &option_index);

     // if the index is -1 then we have reached the end of the options list
     // and we break out of the while loop
     if (c == -1) {
       break;
     }

     // handle the options
     switch (c) {
       case 0:
         break;
       case 'i':
         indexdir = optarg;
         break;
       case 'o':
         outfile = optarg;
         break;
       case 'x':
         x_coord = atoi(optarg);
         break;
       case 'y':
         y_coord = atoi(optarg);
         break;
       case 's':
         score = atof(optarg);
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

  // Make sure an indexes directory is provided
  if (!indexdir) {
    fprintf(stderr, "Please provide the KINC indexes directory (--indexes option).\n");
    exit(-1);
  }
  // Make sure the output directory exists.
  struct stat st = {0};
  if (stat(indexdir, &st) == -1) {
    fprintf(stderr, "The specified indexes directory ,'%s', is missing. Please check the value of the --indexes argument.\n", indexdir);
    exit(-1);
  }

  if (!outfile) {
    fprintf(stderr, "Please provide the output file name (--outfile option).\n");
    exit(-1);
  }

  if (x_coord < 0) {
    fprintf(stderr, "Please provide a positive integer for the x genes (--x option).\n");
    exit(-1);
  }

  if (y_coord < 0) {
    fprintf(stderr, "Please provide a positive integer for the x genes (--x option).\n");
    exit(-1);
  }

  if (score < 0 || score > 1) {
    fprintf(stderr, "Please provide a score between 0 and 1 (--s option).\n");
    exit(-1);
  }


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
  SQLiteQuery query(indexdir);
  query.run(outfile, x_coord, y_coord, score);
}
