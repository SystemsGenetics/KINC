#include "extract.h"


/**
 * Prints the command-line usage instructions for the similarity command
 */
void print_extract_usage() {
  printf("\n");
  printf("Usage: ./kinc extract [options]\n");
  printf("The list of required options:\n");
  printf("  --ematrix|-e     The file name that contains the expression matrix.\n");
  printf("                   The rows must be genes or probe sets and columns are samples\n");
  printf("  --th|-t          The threshold to cut the similarity matrix. Network files will be generated.\n");
  printf("  --method|-m      The correlation method to use. Supported methods include\n");
  printf("                   Pearson's correlation ('pc'), Spearman's rank correlation ('sc')\n");
  printf("                   and Mutual Information ('mi'). Provide either 'pc', 'sc', or\n");
  printf("                   'mi' as values respectively.\n");
  printf("Optional:\n");
  printf("  --headers        Provide this flag if the first line of the ematrix contains\n");
  printf("                   headers.\n");
  printf("  -x               Extract a single similarity value: the x coordinate. Must also use -y\n");
  printf("  -y               Extract a single similarity value: the y coordinate. Must also use -x.\n");
  printf("  --gene1|-1       Extract a single similarity value: The name of the first gene in a singe\n");
  printf("                   pair-wise comparision.  Must be used with --gene2 option.\n");
  printf("  --gene1|-1       Extract a single similarity value: The name of the second gene in a singe\n");
  printf("                   pair-wise comparision.  Must be used with --gene1 option.\n");
  printf("  --max_missing|-g The total number of allowed missing values.  Each gene\n");
  printf("                   comparision can potentially have a missing value in any\n");
  printf("                   sample.  When the maximum number of missing values exceeds\n");
  printf("                   this value the clusters are ignored.\n");
  printf("                   If not provided, no limit is set.\n");
  printf("  --min_csize|-z   The minimum cluster size (number of samples per cluster).\n");
  printf("                   Default is 30\n");
  printf("  --help|-h        Print these usage instructions\n");
  printf("\n");
}

/**
 * The function for extracting results from the similarity matrix
 */
int do_extract(int argc, char *argv[]) {

  // Get the similarity matrix.
//  SimMatrixBinary * smatrix = new SimMatrixBinary(argc, argv);
  SimMatrixTabCluster * smatrix = new SimMatrixTabCluster(argc, argv);

  // If we have a threshold then we want to get the edges of the network.
  // Otherwise the user has asked to print out the similarty value for
  // two genes.
  if (smatrix->getThreshold() > 0) {
    smatrix->writeNetwork();
  }
  else {
    smatrix->getSimilarity();
  }

  return 1;
}

/**
 * Constructor.
 */
SimilarityMatrix::SimilarityMatrix(int argc, char *argv[]) {

  // Set some default values;
  headers = 0;
  rows = 0;
  cols = 0;
  max_missing = INFINITY;
  min_cluster_size = 30;
  strcpy(method, "pc");
  x_coord = -1;
  y_coord = -1;
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
     {"headers",     no_argument,       &headers,  1 },
     {"quiet",       no_argument,       &quiet,  1 },
     {"ematrix",     required_argument, 0,  'e' },
     {"rows",        required_argument, 0,  'r' },
     {"cols",        required_argument, 0,  'c' },
     {"method",      required_argument, 0,  'm' },
     {"th",          required_argument, 0,  't' },
     {"gene1",       required_argument, 0,  '1' },
     {"gene2",       required_argument, 0,  '2' },
     {"max_missing", required_argument, 0,  'g' },
     {"min_csize",   required_argument, 0,  'z' },
     {"help",        no_argument,       0,  'h' },
     {0, 0, 0, 0}  // last element required to be all zeros
   };

   // get the next option
   c = getopt_long(argc, argv, "e:m:x:y:t:h", long_options, &option_index);

   // if the index is -1 then we have reached the end of the options list
   // and we break out of the while loop
   if (c == -1) {
     break;
   }

   // handle the options
   switch (c) {
     case 0:
       break;
     case 'r':
       rows = atoi(optarg);
       break;
     case 'c':
       cols = atoi(optarg);
       break;
     case 'e':
       infilename = optarg;
       break;
     case 'x':
       x_coord = atoi(optarg);
       break;
     case 'y':
       y_coord = atoi(optarg);
       break;
     case 'm':
       strcpy(method, optarg);
       break;
     case 't':
       th = atof(optarg);
       break;
     case '1':
       gene1 = optarg;
       break;
     case '2':
       gene2 = optarg;
       break;
     case 'z':
       min_cluster_size = atoi(optarg);
       break;
     case 'g':
       max_missing = atoi(optarg);
       break;
     case 'h':
       print_extract_usage();
       exit(-1);
       break;
     case '?':
       exit(-1);
       break;
     case ':':
       print_extract_usage();
       exit(-1);
       break;
     default:
       print_extract_usage();
    }
  }

  // make sure the required arguments are set and appropriate
  if (!infilename) {
    fprintf(stderr, "Please provide an expression matrix (--ematrix option). Use the -h option for help.\n");
    exit(-1);
  }

  if (!method) {
    fprintf(stderr, "Please provide the method (--method option) used to construct the similarity matrix.\n");
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

  if (th > 0 && (x_coord > 0 ||  y_coord > 0)) {
    fprintf(stderr, "Please provide a threshold or x and y coordinates only but not both\n");
    exit(-1);
  }

  if (max_missing < -1) {
    fprintf(stderr, "Please provide a positive integer value for maximum missing values\n");
    fprintf(stderr, "the expression matrix (--max_missing option).\n");
    exit(-1);
  }
  if (min_cluster_size < 0 || min_cluster_size == 0) {
    fprintf(stderr, "Please provide a positive integer value for the minimum cluster size in\n");
    fprintf(stderr, "the expression matrix (--min_csize option).\n");
    exit(-1);
  }

  // make sure the input file exists
  if (access(infilename, F_OK) == -1) {
    fprintf(stderr, "The input file does not exists or is not readable.\n");
    exit(-1);
  }

  if (strcmp(method, "pc") != 0 &&
      strcmp(method, "sc") != 0 &&
      strcmp(method, "mi") != 0) {
    fprintf(stderr, "The method (--method option) must either be 'pc' or 'mi'.\n");
    exit(-1);
  }

  // print out some setup details
  if (!quiet) {
    printf("  Using method: '%s'\n", method);
    if (th > 0) {
      printf("  Using threshold of %f\n", th);
    }
  }

  if (strcmp(method, "mi") == 0) {
    input_dir = (char *) malloc(sizeof(char) * 3);
    strcpy(input_dir, "MI");
  }
  else if (strcmp(method, "pc") == 0) {
    input_dir = (char *) malloc(sizeof(char) * 8);
    strcpy(input_dir, "Pearson");
  }
  else if (strcmp(method, "sc") == 0) {
    input_dir = (char *) malloc(sizeof(char) * 9);
    strcpy(input_dir, "Spearman");
  }

  if ((gene1 && !gene2) || (!gene1 && gene2)) {
    fprintf(stderr, "You must provide both gene1 and gene2 options.\n");
    exit(-1);
  }

  // Load the input expression matrix.
  ematrix = new EMatrix(infilename, rows, cols, headers);

  // if the user supplied gene
  if (gene1 && gene2) {
    // Make sure the coordinates are positive integers
    if (x_coord < 0) {
      fprintf(stderr, "Could not find gene %s in the genes list file\n", gene1);
      exit(-1);
    }
    if (y_coord < 0) {
      fprintf(stderr, "Could not find gene %s in the genes list file\n", gene2);
      exit(-1);
    }

    // If the user provided gene names then map those to the numeric IDs.
    char ** genes = ematrix->getGenes();
    int num_genes = ematrix->getNumGenes();
    int i = 0;
    for (i = 0; i < num_genes; i++) {
      if (strcmp(genes[i], gene1) == 0) {
        x_coord = i + 1;
      }
      if (strcmp(genes[i], gene2) == 0) {
         y_coord = i + 1;
      }
    }
  }

  // Make sure we have a positive integer for the x and y coordinates.
  if ((x_coord >= 0 &&  y_coord < 0) ||
      (x_coord < 0  &&  y_coord >= 0)) {
    fprintf(stderr, "Please provide a positive integer for both the x and y coordinates (-x and -y options)\n");
    exit(-1);
  }
}

/**
 * Destructor.
 */
SimilarityMatrix::~SimilarityMatrix() {
  delete ematrix;
}

/**
 * Constructor
 */
SimMatrixBinary::SimMatrixBinary(int argc, char *argv[])
  : SimilarityMatrix(argc, argv){

  // Intialize some variables.
  num_files = 0;

  // Intiialize the file pointer array.
  for(int i = 0; i < 50; i++) {
    files[i] = NULL;
  }

  // Open file handles to all of the binary files.
  openBinFiles();
}

/**
 * Destructor
 */
SimMatrixBinary::~SimMatrixBinary(){

  // Close all of the binary files.
  closeBinFiles();
}

/**
 * Looks in the directory where the similarity matrix is kept and opens each
 * of the bin files for reading.
 */
void SimMatrixBinary::openBinFiles() {
  // A directory handle.
  DIR * FD;
  // A file structure for the file while looping.
  struct dirent* curr_file;
  // the position of the '.bin' extension in the file name.
  char * bin_pos;
  // The position of the .mi or .pc in the file name.
  char * method_pos;
  // The string representation of the number of the bin file.
  char bin_num[5];
  // The numerical bin number.
  int bin_i;
  // Counts the number of files found.
  int num_files = 0;
  // The number of genes in the file.
  int num_genes;

  // Scanning for the files of the similarity matrix
  if (NULL == (FD = opendir(input_dir))) {
    fprintf(stderr, "Error : Failed to open input directory: '%s'\n", input_dir);
    exit(-1);
  }
  while ((curr_file = readdir(FD))) {
    // skip the . and .. dirs
    if (strcmp(curr_file->d_name, ".") == 0 || strcmp(curr_file->d_name, "..") == 0) {
      continue;
    }

    // Make sure this file has a .bin extension.
    bin_pos = strstr(curr_file->d_name, ".bin");
    if(bin_pos) {

      // Make sure the file has the prefix.
      if(strstr(curr_file->d_name, ematrix->getFilePrefix()) != NULL) {
        // Make sure that this file has the method name preceeded by a period.
        char method[5];
        sprintf(method, ".%s", method);
        method_pos = strstr(curr_file->d_name, method);

        if (method_pos) {
          num_files++;

          // Get the numerical value of this bin file.
          int size = (bin_pos - (method_pos + 3)) * sizeof(char);
          memcpy(bin_num, method_pos + (3 * sizeof(char)), size);
          bin_num[size] = 0; // add a terminator to the size
          bin_i = atoi(bin_num);

          // Open the file and store the file handle for later use and store
          // it in the files array using the bin_num as an index.
          char filename[1024];
          sprintf(filename, "%s/%s", input_dir, curr_file->d_name);
          if (!quiet) {
            printf("  Found file: %s\n", filename);
          }
          FILE * fh = fopen(filename, "rb");
          if (fh == NULL) {
            fprintf(stderr, "ERROR: could not open bin file: '%s':\n", filename);
            exit(-1);
          }
          files[bin_i] = fh;
          num_files++;

          // Read in the number of genes in the similarity matrix and the number
          // of lines in the file.
          fread(&num_genes, sizeof(int), 1, fh);
          fread(&num_lines[bin_i], sizeof(int), 1, fh);
        }
      }
    }
  }
  if (num_files == 0) {
    fprintf(stderr, "ERROR: Could not find any matrix .bin files.\n");
    exit(-1);
  }
}

/**
 * Closes all of the file pointers opened by the open_bin_files
 */
void SimMatrixBinary::closeBinFiles() {
  int i;

  for (i = 0; i < num_files; i++) {
    fclose(files[i]);
  }
}

/**
 *
 */
void SimMatrixBinary::writeNetwork() {
   int bin_i = 0;
   int x,y;   // used for iterating through the n x n similarity matrix
   int i;     // used for specifying the bin file to read
   float n;   // the cell value in the similarity matrix

   // the three network output files for edges, negative correlated edges
   // and positive correlated edges
   FILE * edges;
   FILE * edgesN;
   FILE * edgesP;

   // stores the actual name of the output files
   char edges_file[1024];
   char edgesN_file[1024];
   char edgesP_file[1024];

   char * file_prefix = ematrix->getFilePrefix();
   int num_genes = ematrix->getNumGenes();
   char ** genes = ematrix->getGenes();

   sprintf(edges_file, "%s.%s.th%0.6f.coexpnet.edges.txt", file_prefix, method, th);
   edges = fopen(edges_file, "w");
   if (!quiet) {
     printf("  Creating network files...\n");
   }
   fprintf(edges, "gene1\tgene2\tsimilarity\tinteraction\n");


   // The Spearman and Pearson correlation methods will have both negative and
   // positive values, so we want to create separate files for each one.
   if (strcmp(method, "pc") == 0 ||
       strcmp(method, "sc") == 0) {
     sprintf(edgesN_file, "%s.%s.th%0.6f.neg.coexpnet.edges.txt", file_prefix, method, th);
     sprintf(edgesP_file, "%s.%s.th%0.6f.pos.coexpnet.edges.txt", file_prefix, method, th);
     edgesN = fopen(edgesN_file, "w");
     edgesP = fopen(edgesP_file, "w");
     fprintf(edgesN, "gene1\tgene2\tsimilarity\tinteraction\n");
     fprintf(edgesP, "gene1\tgene2\tsimilarity\tinteraction\n");
   }




   // get the size of the matrix in one dimension (i.e. mxm)
   i = 0;
   for (x = 0; x < num_genes; x++) {
      if (i >= num_lines[bin_i]) {
         bin_i++;
         i = 0;
         // set the file position just past the row info
         fseek(files[bin_i], sizeof(int) * 2, SEEK_SET);
      }
      for (y = 0; y < num_genes; y++) {

         // get the next float value for coordinate (x,y)
         int num_read = fread(&n, sizeof(float), 1, files[bin_i]);
         if(num_read != 1){
            fprintf(stderr, "ERROR: cannot fetch from bin file %d\n", bin_i);
            exit(-1);
         }

         // the matrix is symetrical so we don't need to look where y >= x
         if(y >= x) {
           break;
         }

         // write the simiarity value to the appopriate file
         if((n > 0 &&  n >= th) || (n < 0 && -n >= th)){
            fprintf(edges,"%s\t%s\t%0.8f\tco\n", genes[x], genes[y], n);

            // if the method id 'pc' (Pearson's correlation) then we will have
            // negative and positive values, and we'll write those to separate files
            if (strcmp(method, "pc") == 0 ||
                strcmp(method, "sc") == 0) {
              if(n >= 0){
                 fprintf(edgesP, "%s\t%s\t%0.8f\tco\n", genes[x], genes[y], n);
              }
              else {
                fprintf(edgesN, "%s\t%s\t%0.8f\tco\n", genes[x], genes[y], n);
              }
            }
         }
      }
      i++;
   }
   fclose(edges);
   fclose(edgesN);
   fclose(edgesP);
}
/**
 *
 */
void SimMatrixBinary::getPosition() {

  float n;   // the cell value in the similarity matrix
  int temp;  // a temp value used to flip the x & y coordinates if necessary
  int x = x_coord - 1;  // the x coordinate
  int y = y_coord - 1;  // the y coordinate
  int i = 0; // holds the number of lines already visited
  int j;     // iterates through the columns
  int pos;   // used to store the position in the file for the (x,y) coordinate
  int bin_i; // indicates the bin file number currently being checked

  // if y > x then reverse the two as the sim matrix is symetrical and only
  // half is stored.

  if(y > x){
    temp = x;
    x = y;
    y = temp;
  }

  // iterate through the bin files until we find the one with the proper coordinates
  pos = 0;
  for (bin_i = 0; bin_i < num_files; bin_i++) {
    // if the x coordinate falls within the rows stored in the bin file, then
    // find the coordinate value
    if (x >= i && x < i + num_lines[bin_i]) {

      // calculate the coordinate value
      for (j = 0; j < x - i; j++) {
         pos = pos + i + 1 + j;
      }
      pos = (pos + y) * sizeof(float) + sizeof(int) * 2;

      // set the file position to the calculated location of the (x,y) coordinate
      fseek(files[bin_i], pos, SEEK_SET);
      int num_bytes = fread(&n, sizeof(float), 1, files[bin_i]);
      if(num_bytes != 1){
        fprintf(stderr,"ERROR: cannot fetch from bin file %d\n", bin_i);
        exit(-1);
      }
      if (!quiet) {
        printf("similarity(%i,%i) = %0.8f, bin = %d, pos = %d\n", x + 1, y + 1, n, bin_i, pos);
      }
      else {
        printf("%0.8f\n", n);
      }
      break;
    }
    i += num_lines[bin_i];
  }
}

/**
 * Constructor
 */
SimMatrixTabCluster::SimMatrixTabCluster(int argc, char *argv[])
  : SimilarityMatrix(argc, argv) {

  // Initialize the class members.
  num_jobs = 0;

  getNumJobs();
}
/**
 * Destructor
 */
SimMatrixTabCluster::~SimMatrixTabCluster() {

}
/**
 *
 */
void SimMatrixTabCluster::getNumJobs() {
  // Holds the clusters directory name.
  char dirname[1024];

  // Make sure the output directory exists.
  struct stat st = {0};
  char clusterdir[100];
  sprintf(clusterdir, "clusters-%s", method);
  if (stat(clusterdir, &st) == -1) {
    fprintf(stderr, "The clusters directory is missing. Cannot continue.\n");
    exit(-1);
  }

  // Iterate through the files in the 'nan' directory to count the number of
  // jobs used to generate the clusters, we will then intitialize that many
  // file pointers for each of the 102 file arrays.
//  sprintf(dirname, "./%s/%s", clusterdir, "nan");
  sprintf(dirname, "./%s/%d", clusterdir, 100);
  DIR * dir;
  dir = opendir(dirname);
  if (!dir) {
    fprintf(stderr, "The clusters sub directory, %s, is missing. Cannot continue.\n", dirname);
    exit(-1);
  }
  struct dirent * entry;
  while ((entry = readdir(dir)) != NULL) {
    const char * filename = entry->d_name;
    // Skip the . and .. files.
    if (strcmp(filename, ".") == 0 || strcmp(filename, "..") == 0) {
      continue;
    }
    num_jobs++;
  }
}



/**
 *
 */
void SimMatrixTabCluster::writeNetwork() {
   // the three network output files for edges, negative correlated edges
   // and positive correlated edges
   FILE * edges;
   FILE * edgesN;
   FILE * edgesP;

   // stores the actual name of the output files
   char edges_file[1024];
   char edgesN_file[1024];
   char edgesP_file[1024];

   // Get some information from the EMatrix object.
   char * file_prefix = ematrix->getFilePrefix();
   char ** genes = ematrix->getGenes();
   int num_samples = ematrix->getNumSamples();

   // The 8 fields of the clusters tab file.
   int x, y, cluster_num, num_clusters, cluster_samples, num_missing;
   char * samples = (char *) malloc(sizeof(char) * num_samples);
   float cv;

   if (!quiet) {
     printf("  Creating network files...\n");
   }

   // Open the edges output file and write the headers.
   char outfile_prefix[2048];
   if (max_missing > num_samples) {
     sprintf(outfile_prefix, "%s.%s.th%0.6f.mcs%d.mmINF", file_prefix, method, th, min_cluster_size);
   }
   else {
     sprintf(outfile_prefix, "%s.%s.th%0.6f.mcs%d.mm%d", file_prefix, method, th, min_cluster_size, max_missing);
   }
   sprintf(edges_file, "%s.coexpnet.edges.txt", outfile_prefix);
   edges = fopen(edges_file, "w");
   fprintf(edges, "gene1\tgene2\tsimilarity\tinteraction\tcluster\tnum_clusters\tcluster_samples\tmissing_samples\tsamples\n");

   // The Spearman and Pearson correlation methods will have both negative and
   // positive values, so we want to create separate files for each one.
   if (strcmp(method, "pc") == 0 || strcmp(method, "sc") == 0) {
     sprintf(edgesN_file, "%s.neg.coexpnet.edges.txt", outfile_prefix);
     sprintf(edgesP_file, "%s.pos.coexpnet.edges.txt", outfile_prefix);
     edgesN = fopen(edgesN_file, "w");
     edgesP = fopen(edgesP_file, "w");
     fprintf(edgesN, "gene1\tgene2\tsimilarity\tinteraction\tcluster\tnum_clusters\tcluster_samples\tmissing_samples\tsamples\n");
     fprintf(edgesP, "gene1\tgene2\tsimilarity\tinteraction\tcluster\tnum_clusters\tcluster_samples\tmissing_samples\tsamples\n");
   }

   char clusterdir[100];
   char dirname[1024];
   sprintf(clusterdir, "clusters-%s", method);

   // Iterate in descending order through the files and
   int limit = (int) (th * 100.0);
   for (int i = 100; i >= limit; i--) {
     sprintf(dirname, "./%s/%03d", clusterdir, i);
     for (int j = 0; j < num_jobs; j++) {

       // Open the file
       char path[1024];
       sprintf(path, "%s/%s.clusters.%03d.%03d.txt", dirname, ematrix->getFilePrefix(), i, j + 1);
       FILE * fp = fopen(path, "r");
       if (!fp) {
         fprintf(stderr, "Could not open clustering file: '%s'.\n", path);
         exit(-1);
       }

       // Get the values from the line in the file.
       int matches = fscanf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s\n", &x, &y, &cluster_num, &num_clusters, &cluster_samples, &num_missing, &cv, samples);
       while (matches == 8) {
         if (fabs(cv) >= th && cluster_samples >= min_cluster_size  && num_missing <= max_missing) {
           fprintf(edges, "%s\t%s\t%0.8f\tco\t%d\t%d\t%d\t%d\t%s\n", genes[x-1], genes[y-1], cv, cluster_num, num_clusters, cluster_samples, num_missing, samples);

           // if the method id 'pc' (Pearson's correlation) then we will have
           // negative and positive values, and we'll write those to separate files
           if (strcmp(method, "pc") == 0 || strcmp(method, "sc") == 0) {
             if(cv >= 0){
                fprintf(edgesP, "%s\t%s\t%0.8f\tco\t%d\t%d\t%d\t%d\t%s\n", genes[x-1], genes[y-1], cv, cluster_num, num_clusters, cluster_samples, num_missing, samples);
             }
             else {
               fprintf(edgesN, "%s\t%s\t%0.8f\tco\t%d\t%d\t%d\t%d\t%s\n", genes[x-1], genes[y-1], cv, cluster_num, num_clusters, cluster_samples, num_missing, samples);
             }
           }
         }
         matches = fscanf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s\n", &x, &y, &cluster_num, &num_clusters, &cluster_samples, &num_missing, &cv, samples);
       }
       fclose(fp);
     }
   }

   free(samples);
   fclose(edges);
   fclose(edgesN);
   fclose(edgesP);
}
