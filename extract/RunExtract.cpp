#include "RunExtract.h"

/**
 * Prints the command-line usage instructions for the similarity command
 */
void RunExtract::printUsage() {
  printf("\n");
  printf("Usage: ./kinc extract [options]\n");
  printf("The list of required options:\n");
  printf("  --ematrix|-e     The file name that contains the expression matrix.\n");
  printf("                   The rows must be genes or probe sets and columns are samples\n");
  printf("  --rows|-r        The number of lines in the ematrix file including the header\n");
  printf("                   row if it exists\n");
  printf("  --cols|-c        The number of columns in the input file\n");
  printf("  --th|-t          The threshold to cut the similarity matrix. Network files will be generated.\n");
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
  printf("Optional filtering arguments:\n");
  printf("  -x               Extract a single similarity value: the x coordinate. Must also use -y\n");
  printf("  -y               Extract a single similarity value: the y coordinate. Must also use -x.\n");
  printf("  --gene1|-1       Extract a single similarity value: The name of the first gene in a singe\n");
  printf("                   pair-wise comparision.  Must be used with --gene2 option.\n");
  printf("  --gene1|-1       Extract a single similarity value: The name of the second gene in a singe\n");
  printf("                   pair-wise comparision.  Must be used with --gene1 option.\n");
  printf("\n");
  printf("Optional arguments for clustered data:\n");
  printf("  --clustering|-l  The type of clustering that was performed during construction\n");
  printf("                   of the similarity matrix (e.g. 'mixmod').\n");
  printf("  --max_missing|-g The total number of allowed missing values.  Each gene\n");
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
  printf("  --help|-h        Print these usage instructions\n");
  printf("\n");
}

RunExtract::RunExtract(int argc, char *argv[]) {
  // Set some default values;
   headers = 0;
   rows = 0;
   cols = 0;
   max_missing = INFINITY;
   max_modes = 1;
   min_cluster_size = 30;
   method = NULL;
   x_coord = -1;
   y_coord = -1;
   gene1 = NULL;
   gene2 = NULL;
   th = 0;
   quiet = 0;

   // The value returned by getopt_long
   int c;

   // loop through the incoming arguments until the
   // getopt_long function returns -1. Then we break out of the loop
   while(1) {
    int option_index = 0;

    // specify the long options. The values returned are specified to be the
    // short options which are then handled by the case statement below
    static struct option long_options[] = {
      {"quiet",        no_argument,       &quiet,  1 },
      {"method",       required_argument, 0,  'm' },
      {"help",         no_argument,       0,  'h' },
      // Expression matrix options.
      {"rows",         required_argument, 0,  'r' },
      {"cols",         required_argument, 0,  'c' },
      {"headers",      no_argument,       &headers,  1 },
      {"omit_na",      no_argument,       &omit_na,  1 },
      {"func",         required_argument, 0,  'f' },
      {"na_val",       required_argument, 0,  'n' },
      {"ematrix",      required_argument, 0,  'e' },
      // Common fitering options
      {"th",           required_argument, 0,  't' },
      {"gene1",        required_argument, 0,  '1' },
      {"gene2",        required_argument, 0,  '2' },
      {"x",            required_argument, 0,  'x' },
      {"y",            required_argument, 0,  'y' },
      // Clustered options
      {"clustering",   required_argument, 0,  'l' },
      {"max_missing",  required_argument, 0,  'g' },
      {"min_csize",    required_argument, 0,  'z' },
      {"max_modes",    required_argument, 0,  'd' },

      // Last element required to be all zeros.
      {0, 0, 0, 0}
    };
    delete ematrix;

    // get the next option
    c = getopt_long(argc, argv, "m:r:c:f:n:e:t:1:2:x:y:g:d:z:l:h", long_options, &option_index);

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
      // Clustering options.
      case 'l':
        clustering = optarg;
        break;
      // Common fitering options
      case 't':
        th = atof(optarg);
        break;
      case '1':
        gene1 = optarg;
        break;
      case '2':
        gene2 = optarg;
        break;
      case 'x':
        x_coord = atoi(optarg);
        break;
      case 'y':
        y_coord = atoi(optarg);
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
      // Clustered data filters.
      case 'z':
        min_cluster_size = atoi(optarg);
        break;
      case 'g':
        max_missing = atoi(optarg);
        break;
      case 'd':
        max_modes = atoi(optarg);
        break;
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

   if (th > 0 && (x_coord > 0 ||  y_coord > 0)) {
     fprintf(stderr, "Please provide a threshold or x and y coordinates only but not both\n");
     exit(-1);
   }

   if (max_missing < -1) {
     fprintf(stderr, "Please provide a positive integer value for maximum missing values\n");
     fprintf(stderr, "the expression matrix (--max_missing option).\n");
     exit(-1);
   }
   if (max_modes < 1) {
     fprintf(stderr, "Please provide a positive integer greater than 1 for the maximum modes values (--max_modes option).\n");
     exit(-1);
   }
   if (min_cluster_size < 0 || min_cluster_size == 0) {
     fprintf(stderr, "Please provide a positive integer value for the minimum cluster size in\n");
     fprintf(stderr, "the expression matrix (--min_csize option).\n");
     exit(-1);
   }

   if ((gene1 && !gene2) || (!gene1 && gene2)) {
     fprintf(stderr, "You must provide both gene1 and gene2 options.\n");
     exit(-1);
   }

   // Load the input expression matrix.
   ematrix = new EMatrix(infilename, rows, cols, headers, omit_na, na_val, func);

   // if the user supplied gene
   if (gene1 && gene2) {
     x_coord = findGeneCoord(gene1);
     y_coord = findGeneCoord(gene2);

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

   // print out some setup details
   if (!quiet) {
     printf("  Using method: '%s'\n", method);
     if (th > 0) {
       printf("  Using threshold of %f\n", th);
     }
     else {
       printf("  Using coords (%d, %d)\n", x_coord, y_coord);
     }
   }
}
/**
 *
 */
RunExtract::~RunExtract() {
  delete ematrix;
}

/**
 *
 */
int RunExtract::findGeneCoord(char * gene) {
  char ** genes = ematrix->getGenes();
  int num_genes = ematrix->getNumGenes();

  for (int i = 1; i <= num_genes; i++) {
    if (strcmp(gene, genes[i-1]) == 0) {
      return i;
    }
  }
  return -1;
}

/**
 *
 */
void RunExtract::execute() {

  // Get the similarity matrix.
  if (clustering) {
    SimMatrixTabCluster * smatrix = new SimMatrixTabCluster(ematrix, quiet,
      method, x_coord, y_coord, gene1, gene2, th, max_missing, min_cluster_size,
      max_modes);
    // If we have a threshold then we want to get the edges of the network.
    // Otherwise the user has asked to print out the similarity value for
    // two genes.
    if (smatrix->getThreshold() > 0) {
      smatrix->writeNetwork();
    }
    else {
      smatrix->getPosition();
    }
  }
  else {
    SimMatrixBinary * smatrix = new SimMatrixBinary(ematrix, quiet, method,
      x_coord, y_coord, gene1, gene2, th);
    // If we have a threshold then we want to get the edges of the network.
    // Otherwise the user has asked to print out the similarity value for
    // two genes.
    if (smatrix->getThreshold() > 0) {
      smatrix->writeNetwork();
    }
    else {
      smatrix->getPosition();
    }
  }

}
