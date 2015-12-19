#include "SimMatrixBinary.h"

/**
 * Constructor
 */
SimMatrixBinary::SimMatrixBinary(EMatrix *ematrix, int quiet, char ** method, int num_methods,
    char * th_method, int x_coord, int y_coord, char * gene1, char * gene2, float th)
  : SimilarityMatrix(ematrix, quiet, method, num_methods, th_method, x_coord, y_coord, gene1, gene2, th){

  // Initialize some variables.
  num_files = 0;

  // Initialize the file pointer array.
  for(int i = 0; i < 50; i++) {
    files[i] = NULL;
  }

  // For the binary file format:
  bin_dir = (char *) malloc(sizeof(char) * strlen(ematrix->getInfileName()));
  if (strcmp(th_method, "mi") == 0) {
    strcpy(bin_dir, "MI");
  }
  else if (strcmp(th_method, "pc") == 0) {
    strcpy(bin_dir, "Pearson");
  }
  else if (strcmp(th_method, "sc") == 0) {
    strcpy(bin_dir, "Spearman");
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
  free(bin_dir);
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
  // The number of genes in the file.
  int num_genes;

  // Scanning for the files of the similarity matrix
  if (!quiet) {
    printf("  Scanning %s directory for matrix .bin files...\n", bin_dir);
  }
  if (NULL == (FD = opendir(bin_dir))) {
    fprintf(stderr, "Error : Failed to open input directory: '%s'\n", bin_dir);
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

        // Make sure that this file has the method name proceeded by a period.
        char m[5];
        sprintf(m, ".%s", th_method);
        method_pos = strstr(curr_file->d_name, m);

        if (method_pos) {

          // Get the numerical value of this bin file.
          int size = (bin_pos - (method_pos + 3)) * sizeof(char);
          memcpy(bin_num, method_pos + (3 * sizeof(char)), size);
          bin_num[size] = 0; // add a terminator to the size
          bin_i = atoi(bin_num);

          // Open the file and store the file handle for later use and store
          // it in the files array using the bin_num as an index.
          char filename[1024];
          sprintf(filename, "%s/%s", bin_dir, curr_file->d_name);
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

   sprintf(edges_file, "%s.%s.th%0.6f.coexpnet.edges.txt", file_prefix, th_method, th);
   edges = fopen(edges_file, "w");
   if (!quiet) {
     printf("  Creating network files...\n");
   }
   fprintf(edges, "gene1\tgene2\tsimilarity\tinteraction\n");


   // The Spearman and Pearson correlation methods will have both negative and
   // positive values, so we want to create separate files for each one.
   if (strcmp(th_method, "pc") == 0 ||
       strcmp(th_method, "sc") == 0) {
     sprintf(edgesN_file, "%s.%s.th%0.6f.neg.coexpnet.edges.txt", file_prefix, th_method, th);
     sprintf(edgesP_file, "%s.%s.th%0.6f.pos.coexpnet.edges.txt", file_prefix, th_method, th);
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
            if (strcmp(th_method, "pc") == 0 ||
                strcmp(th_method, "sc") == 0) {
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

  // The cell value in the similarity matrix.
  float n;
  // A temporary value used to flip the x & y coordinates if necessary.
  int temp;
  // The x coordinate.
  int x = x_coord - 1;
  // The y coordinate.
  int y = y_coord - 1;
  // Holds the number of lines already visited.
  int i = 0;
  // Iterates through the columns.
  int j;
  // Used to store the position in the file for the (x,y) coordinate.
  int pos;
  // Indicates the bin file number currently being checked.
  int bin_i;

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
