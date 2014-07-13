#include "extract.h"

/**
 * The function for extracting results from the similarity matrix
 */
int do_extract(int argc, char *argv[]) {

  int c;  // The value returned by getopt_long

  static NetParameters params;

  // initialize some of the program parameters
  strcpy(params.method, "pc");
  params.x_coord = 0;
  params.y_coord = 0;
  params.th = 0;
  params.num_files = 0;

  // loop through the incoming arguments until the
  // getopt_long function returns -1. Then we break out of the loop
  while(1) {
    int option_index = 0;

    // specify the long options. The values returned are specified to be the
    // short options which are then handled by the case statement below
    static struct option long_options[] = {
      {"ematrix", required_argument, 0,  'e' },
      {"method",  required_argument, 0,  'm' },
      {"th",      required_argument, 0,  't' },
      {"genes",   required_argument, 0,  'g' },
      {"gene1",   required_argument, 0,  '1' },
      {"gene2",   required_argument, 0,  '2' },
      {0, 0, 0,  0 }  // last element required to be all zeros
    };

    // get the next option
    c = getopt_long(argc, argv, "e:m:x:y:t", long_options, &option_index);

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
      case 'g':
        params.genes_file = optarg;
        break;
      case 'x':
        params.x_coord = atoi(optarg);
        break;
      case 'y':
        params.y_coord = atoi(optarg);
        break;
      case 'm':
        strcpy(params.method, optarg);
        break;
      case 't':
        params.th = atof(optarg);
        break;
      case '1':
        params.gene1 = optarg;
        break;
      case '2':
        params.gene2 = optarg;
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

  // make sure the required arguments are set and appropriate
  if (!params.infilename) {
    fprintf(stderr,"Please provide an expression matrix (--ematrix option). Use the -h option for help.\n");
    exit(-1);
  }

  if (!params.method) {
    fprintf(stderr,"Please provide the method (--method option) used to construct the similarity matrix.\n");
    exit(-1);
  }

  if (!params.genes_file) {
    fprintf(stderr,"Please provide the genes file (--genes option).\n");
    exit(-1);
  }

  // make sure we have a positive integer for the rows and columns of the matrix
  if ((params.x_coord >  0 && params.y_coord == 0) ||
      (params.x_coord == 0 && params.y_coord  > 0)) {
    fprintf(stderr,"Please provide a positive integer for both the x and y coordinates (-x and -y options)\n");
    exit(-1);
  }

  if (params.th > 0 && (params.x_coord > 0 || params.y_coord > 0)) {
    fprintf(stderr,"Please provide a threshold or x and y coordinates only but not both\n");
    exit(-1);
  }

  // make sure the input file exists
  if (access(params.infilename, F_OK) == -1) {
    fprintf(stderr,"The input file does not exists or is not readable.\n");
    exit(-1);
  }

  // make sure the input file exists
  if (access(params.genes_file, F_OK) == -1) {
    fprintf(stderr,"The genes file does not exists or is not readable.\n");
    exit(-1);
  }

  if (strcmp(params.method, "pc") != 0 && strcmp(params.method, "mi") != 0) {
    fprintf(stderr,"The method (--method option) must either be 'pc' or 'mi'.\n");
    exit(-1);
  }

  // print out some setup details
  printf("  Using method: '%s'\n", params.method);
  if (params.th > 0) {
    printf("  Using threshold of %f\n", params.th);
  }

  // remove the path and extension from the filename
  char * temp = basename(params.infilename);
  strcpy(params.fileprefix, temp);
  char * p = rindex(params.fileprefix, '.');
  if (p) {
    p[0] = 0;
  }

  if (strcmp(params.method, "mi") == 0) {
    params.inputDir = "MI";
  }
  else if (strcmp(params.method, "pc") == 0) {
    params.inputDir = "Pearson";
  }

  // open all of the bin files for easy access
  open_bin_files(&params);

  // get the list of gene names
  get_gene_names(params);

  if (params.th > 0) {
    get_edges(params);
  }
  else {
    printf("  Finding similarity value at position (%d, %d)\n", params.x_coord, params.y_coord);
    get_position(params);
  }

  // close all of the bin files
  close_bin_files(params);

  return 1;
}

/**
 * Looks in the directory where the similarity matrix is kept and opens each
 * of the bin files for reading.
 */
void open_bin_files(NetParameters *params) {
  DIR * FD;                  // a directory handle
  struct dirent* curr_file;  // a file structure for the file while looping
  char * bin_pos;            // the position of the '.bin' extension in the file name
  char * method_pos;         // the position of the .mi or .pc in the file name
  char bin_num[5];           // the string representation of the number of the bin file

  // Scanning for the files of the similarity matrix
  if (NULL == (FD = opendir(params->inputDir))) {
    fprintf(stderr, "Error : Failed to open input directory - %s\n", strerror(errno));
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
      if(strstr(curr_file->d_name, params->fileprefix) != NULL) {

        // Make sure that this file has the method name preceeded by a period.
        char method[5];
        sprintf(method, ".%s", params->method);
        method_pos = strstr(curr_file->d_name, method);
        if (method_pos) {
          // Get the numerical value of this bin file.
          int size = (bin_pos - (method_pos + 3)) * sizeof(char);
          memcpy(bin_num, method_pos + (3 * sizeof(char)), size);
          bin_num[size] = 0; // add a terminator to the size

          // Open the file and store the file handle for later use and store
          // it in the files array using the bin_num as an index.
          printf("  Opening file: %s\n", curr_file->d_name);
          params->files[atoi(bin_num)] = fopen(curr_file->d_name, "w");
          params->num_files++;

          // Read in the number of genes in the similarity matrix and the number
          // of lines in the file.
          fread(&params->numGenes, sizeof(int), 1, params->files[atoi(bin_num)]);
          fread(&params->numLines[atoi(bin_num)], sizeof(int), 1, params->files[atoi(bin_num)]);
          printf(" NUM GENES: %d\n", params->numGenes);
          printf(" NUM LINES: %d, %d\n", params->numLines[atoi(bin_num)], atoi(bin_num));
        }
      }
    }
  }
}

/**
 * Closes all of the file pointers opened by the open_bin_files
 */
void close_bin_files(NetParameters params) {
  int i;

  for (i = 0; i < params.num_files; i++) {
    fclose(params.files[i]);
  }
}

/**
 * Read in the gene names of the smatrix
 */
void get_gene_names(NetParameters params) {

   FILE * genesf = fopen(params.genes_file, "r");
   char * gene;
   int i = 0;

   // reserve the memory for the genes array (array of strings)
   params.genes = (char **) malloc(sizeof(char *) * params.numGenes);

   while (fscanf(genesf, "%s", params.genes[i++])) {
//      // if $x and $y are not numeric then we assume they are probeset
//      // names and we need to look them up
//      if($x and $x eq $_){
//         $x = $i;
//      }
//      if($y and $y eq $_){
//         $y = $i;
//      }
   }
   fclose(genesf);
}

/**
 *
 */
void get_edges(NetParameters params) {
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

   sprintf(edges_file, "%s.coexpnet.edges.txt", params.fileprefix);
   sprintf(edgesN_file, "%s.neg.coexpnet.edges.txt", params.fileprefix);
   sprintf(edgesP_file, "%s.pos.coexpnet.edges.txt", params.fileprefix);

   edges = fopen(edges_file, "w");
   edgesN = fopen(edgesN_file, "w");
   edgesP = fopen(edgesP_file, "w");

   // get the size of the matrix in one dimension (i.e. mxm)
   i = 0;
   for (x = 0; x < params.numGenes; x++) {
      if (i >= params.numLines[bin_i]) {
         bin_i++;
         i = 0;
         // set the file position just past the row info
         fseek(params.files[bin_i], 8, SEEK_SET);
      }
      for (y = 0; y < params.numGenes; y++) {

         // get the next float value for coordinate (x,y)
         size_t num_bytes = fread(&n, 4, 1, params.files[bin_i]);
         if(num_bytes != 4){
            fprintf(stderr,"ERROR: cannot fetch from bin file $bin_i\n");
            exit(-1);
         }

         // the matrix is symetrical so we don't need to look where y >= x
         if(y >= x) {
           break;
         }

         // write the simiarity value to the appopriate file
         if((n > 0 &&  n >= params.th) || (n < 0 && -n >= params.th)){
            fprintf(edges,"%s\t%s\t%0.8f\n", params.genes[x], params.genes[y], n);

            // if the method id 'pc' (Pearson's correlation) then we will have
            // negative and positive values, and we'll write those to separate files
            if (strcmp(params.method, "pc")) {
              if(n >= 0){
                 fprintf(edgesP, "%s\t%s\t%0.8f\n", params.genes[x], params.genes[y], n);
              }
              else {
                fprintf(edgesN, "%s\t%s\t%0.8f\n", params.genes[x], params.genes[y], n);
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
void get_position(NetParameters params) {

  float n;   // the cell value in the similarity matrix
  int temp;  // a temp value used to flip the x & y coordinates if necessary
  int x = params.x_coord - 1;  // the x coordinate
  int y = params.y_coord - 1;  // the y coordinate
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


  printf(" num_files: %d\n", params.num_files);
  for (bin_i = 0; bin_i < params.num_files; bin_i++) {
    // if the x coordinate falls within the rows stored in the bin file, then
    // find the coordinate value
    printf("x: %d, i: %d, %d\n", x, i, params.numLines[bin_i]);
    if (x >= i && x < i + params.numLines[bin_i]) {


      // calculate the coordinate value
      for (j = 0; j < x - i; j++) {
         pos = pos + i + 1 + j;
      }
      pos = (pos + y) * 4 + 8;

      // set the file position to the calculated location of the (x,y) coordinate
      fseek(params.files[bin_i], pos, SEEK_SET);
      size_t num_bytes = fread(&n, 4, 1, params.files[bin_i]);
      if(num_bytes != 4){
        fprintf(stderr,"ERROR: cannot fetch from bin file $bin_i\n");
        exit(-1);
      }
      printf("cor(%i,%i) = %0.8f, bin = $bin_i, pos = %d\n", x + 1, y + 1, n, pos);
      break;
    }
    i += params.numLines[bin_i];
  }
}

/**
 * Prints the command-line usage instructions for the similarity command
 */
void print_extract_usage() {
  printf("\n");
  printf("Usage: ./RMTGeneNet extract [options]\n");
  printf("The list of required options:\n");
  printf("  --ematrix|-e The file name that contains the expression matrix.\n");
  printf("                 The rows must be genes or probe sets and columns are samples\n");
  printf("  --th|-t      The threshold to cut the similarity matrix. Network files will be generated.\n");
  printf("  --method|-m  The correlation method to use. Supported methods include\n");
  printf("                 Pearson's correlation and Mutual Information. Provide\n");
  printf("                 either 'pc' or mi' as values respectively.\n");
  printf("  --genes|-g   The file containing the list of genes, in order, from the similarity matrix\n");
  printf("  -x           Extract a single similarity value: the x coordinate. Must also use -y\n");
  printf("  -y           Extract a single similarity value: the y coordinate. Must also use -x.\n");
  printf("  --gene1|-1   Extract a single similarity value: The name of the first gene in a singe\n");
  printf("                 pair-wise comparision.  Must be used with --gene2 option.\n");
  printf("  --gene1|-1   Extract a single similarity value: The name of the second gene in a singe\n");
  printf("                 pair-wise comparision.  Must be used with --gene1 option.\n");
  printf("\n");
}
