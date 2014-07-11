#include "threshold.h"

/**
 * The function to call when running the 'similarity' command.
 */
int do_threshold(int argc, char *argv[]) {

  time_t start_time, end_time;  // variables used for timing of the software
  int c;                        // the value returned by getopt_long

  static RMTParameters params;

  // initialize some of the program parameters
  params.perf = 1;
  params.thresholdStart = 0.99;
  params.thresholdStep  = 0.001;
  params.chiSoughtValue = 200;
  params.min_size = 100;
  strcpy(params.method, "pc");

  // loop through the incoming arguments until the
  // getopt_long function returns -1. Then we break out of the loop
  while(1) {
    int option_index = 0;

    // specify the long options. The values returned are specified to be the
    // short options which are then handled by the case statement below
    static struct option long_options[] = {
      {"perf",    no_argument,       &params.perf,     1 },
      {"ematrix", required_argument, 0,  'e' },
      {"method",  required_argument, 0,  'm' },
      {"th",      required_argument, 0,  't' },
      {"chi",     required_argument, 0,  'c' },
      {"step",    required_argument, 0,  's' },
      {0, 0, 0,  0 }  // last element required to be all zeros
    };

    // get the next option
    c = getopt_long(argc, argv, "e:m:t:c:s:", long_options, &option_index);

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
      case 'm':
        strcpy(params.method, optarg);
        printf("  Using method: '%s'\n", params.method);
        break;
      case 't':
        params.thresholdStart = atof(optarg);
        printf("  Start treshold: %f\n", params.thresholdStart);
        break;
      case 'c':
        params.chiSoughtValue = atof(optarg);
        printf("  Desired Chi-square %f\n", params.chiSoughtValue);
        break;
      case 's':
        params.thresholdStep = atof(optarg);
        printf("  Step per iteration: %f\n", params.thresholdStep);
        break;
      case '?':
        exit(-1);
        break;
      case ':':
        print_threshold_usage();
        exit(-1);
        break;
      default:
        print_threshold_usage();
    }
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
  char * temp = basename(params.infilename);
  strcpy(params.fileprefix, temp);
  char * p = rindex(params.fileprefix, '.');
  if (p) {
    p[0] = 0;
  }

  // if performance monitoring is enabled the set the start time
  if (params.perf) {
    time(&start_time);
  }

  if (strcmp(params.method, "mi") == 0) {
    params.inputDir = "MI";
  }
  else if (strcmp(params.method, "pc") == 0) {
    params.inputDir = "Pearson";
  }

  // open the file and get the number of genes and the lines per file
  // these data are the first two integers in the file
  char filename[1024];
  FILE* info;
  sprintf(filename, "%s/%s.%s%d.bin", params.inputDir, params.fileprefix, params.method, 0);
  // TODO: check that file exists before trying to open
  info = fopen(filename, "rb");
  fread(&params.numGenes, sizeof(int), 1, info);
  fread(&params.numLinesPerFile, sizeof(int), 1, info);
  fclose(info);

  printf("  Genes: %d, Num lines per file: %d\n", params.numGenes, params.numLinesPerFile);

  //initialize the global RMTParameters struct
  params.nnsdHistogramBin       = 0.05;
  params.chiSquareTestThreshold = 99.607;
  params.minUnfoldingPace       = 10;
  params.maxUnfoldingPace       = 41;
  params.mimiModuleSize         = 4;
  params.edHistogramBin         = 0.1;

  // allocate memory for s and cutM_index arrays
  params.UsedFlag = (int *) malloc(params.numGenes * sizeof(int));
  params.cutM_index = (int *) malloc(params.numGenes * sizeof(int));

  find_threshold(params);

  // free memory
  free(params.cutM_index);
  free(params.UsedFlag);

  // if performance monitoring is enabled then write the timing data
  if(params.perf) {
    time(&end_time);
    FILE * timingfile = fopen("timingdata.txt", "a");
    fprintf(timingfile, "Minutes in determining: %f\n", (end_time - start_time) / 60.0);
  }

  printf("Done.\n");
  return 1;
}

/**
 * writes a value to a string, then returns the pointer.
 * Only supports base 10 numbers.
 */
char* itoa(int val, char* ptr){
  sprintf(ptr,"%d",val);
  return ptr;
}

/*
 *
 *
 */
int find_threshold(RMTParameters params) {

  float* newM;
  int size;
  time_t start, end;
  FILE* eigenF, *chiF;
  char chi_filename[1024];  // the output file name
  char eigen_filename[1024];  // the output file name

  float th = params.thresholdStart;
  double finalTH  = 0.0;
  double finalChi = 10000.0;
  double minTH    = 1.0;
  double minChi   = 10000.0;
  double maxTH    = 0.0;
  double maxChi   = 0.0;
  int i           = 0;
  double chi;
  float * E;  //array for eigenvalues

  // open the output files
  sprintf(chi_filename, "%s.chiVals.txt", params.fileprefix);
  chiF = fopen(chi_filename, "w");
  sprintf(eigen_filename, "%s.eigenVals.txt", params.fileprefix);
  eigenF = fopen(eigen_filename, "w");
  fprintf(chiF, "Threshold\tChi-square\tCut Matrix Size\n");

  do {
    // decrement the threshold using the step value and then retreive the 
    // matrix that contains only the threshold and higher
    th = th - params.thresholdStep;
    printf("  testing threshold: %f...\n", th);
    printf("  reading bin files...\n");
    newM = read_similarity_matrix(th, &size, params);

    printf("  found matrix of size n x n, n = %d...\n", size);
    if (size >= params.min_size) {
      printf("  calculating eigenvalues...\n");
      E = calculateEigen(newM, size);
      free(newM);

      // print out eigenvalues to file
      fprintf(eigenF, "%f\t", th);
      for (i = 0; i < size ; i++) {
        fprintf(eigenF, "%f\t", E[i]);
      }
      fprintf(eigenF,"\n");

      printf("  testing similarity of NNSD with Poisson...\n");
      chi = chiSquareTestUnfoldingNNSDWithPoisson(E, size, params);
      free(E);

      // if the chi-square test did not fail (== -1) then set the values
      // for the next iteration
      if (chi != -1) {
        fprintf(chiF, "%f\t%f\t%d\n", th, chi, size);
        fflush(chiF);
        printf("  chi = %f\n", chi);

        if(chi < minChi){
          minChi = chi;
          minTH = th;
        }
        if (chi < params.chiSquareTestThreshold){
          finalTH = th;
          finalChi = chi;
        }
        if (finalChi < params.chiSquareTestThreshold && chi > finalChi && th < finalTH){
          maxChi = chi;
          maxTH = th;
        }
      }
    }
    else{
      free(newM);
    }
  }
  while(maxChi < params.chiSoughtValue || size == params.numGenes);


  // If finalChi is still greater than threshold, check the small scale
  if (finalChi > params.chiSquareTestThreshold) {
    fprintf(chiF, "checking small scale\n");
    th = (float)minTH + 0.2;
    for (i = 0 ; i <= 40 ; i++) {
      th = th - params.thresholdStep * i;
      newM = read_similarity_matrix(th, &size, params);

      if (size >= 100) {
        E = calculateEigen(newM, size);
        free(newM);

        // print out eigenvalues to file
        fprintf(eigenF, "%f\t", th);
        for (i=0 ; i<size ; i++) {
          fprintf(eigenF, "%f\t", E[i]);
        }
        fprintf(eigenF, "\n");
        chi = chiSquareTestUnfoldingNNSDWithPoisson(E, size, params);
        fprintf(chiF, "%f\t%f\t%d\n", th, chi, size);
        fflush(chiF);
        free(E);

        if (chi < minChi) {
          minChi = chi;
          minTH = th;
        }
        if (chi < params.chiSquareTestThreshold) {
          finalTH = th;
          finalChi = chi;
        }
      } // end if size >= 100
      else{
        free(newM);
      }
    } // end for 1 -> 40 loop
  } // end if finalChi > rmt...

  // close the chi and eigen files now that results are written
  fclose(chiF);
  fclose(eigenF);

  // Set the Properties file according to success or failure
  if(finalChi < params.chiSquareTestThreshold){
    finalTH = ceil(finalTH * 10000) / 10000.0;
    FILE* th;
    char filename[1024];
    sprintf(filename, "%s.th.txt", params.fileprefix);
    th = fopen(filename, "w");
    fprintf(th, "%f", finalTH);
    fclose(th);
    return 0;
  }
  else{
    finalTH = ceil(finalTH * 10000) / 10000.0;
    return -2;
  }
}

/*
 * Parses the correlation bin files to find pairs of genes with a 
 * correlation value greater than the given theshold and constructs
 * a new correlation matrix that only contains those pairs of genes.
 *
 * @param float th
 *  The minimum threshold to search for. 
 * @param int* size
 *  The size, n, of the cut n x n matrix. This value gets set by the function.
 *
 * @return 
 *  A pointer to a floating point array.  The array is a correlation
 *  matrix containing only the genes that have at least one correlation value 
 *  greater than the given threshold.
 */

float * read_similarity_matrix(float th, int * size, RMTParameters params) {

  float * cutM;    // the resulting cut similarity matrix
  float * rowi;    // holds the float value from row i in the bin file
  int len;         // the length of the filename and path to the correlation bin file
  int i, h;        // used to iterate through the bin files
  int j;           // used to iterate through the rows of each bin file
  int k;           // used to iterate through the cols of each bine file
  int z;           // the number of binary files
  int used;        // holds the number of genes (probesets) that have a greater thrshold
  int limit;       // the maximum row in the current correlation bin file
  int junk;        // a dummy variable
  FILE* in;
  char filename[1024]; // used for storing the bin file name

  memset(params.cutM_index, -1, sizeof(int) * (params.numGenes));

  rowi = (float*) malloc(sizeof(float) * params.numGenes);

  // we need to know how many rows and columns we will have in our cut matrix.
  // the cut matrix is the matrix that only contains genes with a threshold
  // value greater than the given value.  Therefore, we iterate through the 
  // file to find out how many genes we will have in the cut matrix, then
  // we iterate again to build the cut matrix.

  // TODO: we should save an index in the file for where these correlation
  // values are stored rather than iterating through the file twice.

  // ------------------------------
  // Step #1: iterate through the binary correlation files to find out how many
  // genes there will be in the cut matrix.
  z = (params.numGenes - 1) / params.numLinesPerFile;
  for (i = 0; i <= z; i++) {

    sprintf(filename, "%s/%s.%s%d.bin", params.inputDir, params.fileprefix, params.method, i);
    in = fopen(filename, "rb");
    fread(&junk, sizeof(int), 1, in); // numGenes
    fread(&junk, sizeof(int), 1, in); // numLinesPerFile
    if (i != z) {
      limit = (i + 1) * params.numLinesPerFile;
    }
    else{
      limit = params.numGenes;
    }

    // iterate through the rows and columns of the file and look for 
    // entries greater than the provided threshold.  When found, use
    // the row and column indexes to set a '1' in the  array.
    // this array indicates which genes have values we want to keep. 
    for (j = i * params.numLinesPerFile; j < limit; j++) {
      fread(rowi, sizeof(float), j + 1, in);
      for (k = 0; k < j + 1; k++) {
        // if the correlation value is greater than the given threshold then 
        // flag the row/column indexes
        if (k != j && fabs(rowi[k]) > th) {
          params.UsedFlag[k] = 1;
          params.UsedFlag[j] = 1;
        }
      }
    }
    fclose(in);
  }

  // get the number of genes (or probe sets) that have a correlation value
  // greater than the provided threshold value
  used = 0;
  j = 0;
  for (i = 0; i < params.numGenes; i++) {
    if (params.UsedFlag[i] == 1) {
      used++;
      params.cutM_index[i] = j;
      j++;
    }
  }

  // now that we know how many genes have a threshold greater than the
  // given we can allocate memory for new cut matrix
  cutM = (float *) calloc(used * used, sizeof(float));
  // initialize the diagonal to 1
  for (i = 0; i < used; i++) {
    cutM[i + i * used] = 1;
  }

  // set the incoming size argument to be the size dimension of the cut matrix
  *size = used;

  // ------------------------------
  // Step #2: Now build the cut matrix by retrieving the correlation values
  // for each of the genes identified previously.
  for (i = 0; i <= z; i++) {

    sprintf(filename, "%s/%s.%s%d.bin", params.inputDir, params.fileprefix, params.method, i);
    in = fopen(filename, "rb");
    fread(&junk, sizeof(int), 1, in); // numGenes
    fread(&junk, sizeof(int), 1, in); // numLinesPerFile
    if (i != z) {
      limit = (i + 1) * params.numLinesPerFile;
    }
    else{
      limit = params.numGenes;
    }
    // iterate through the rows of the bin file
    for (j = i * params.numLinesPerFile; j < limit; j++) {
      fread(rowi, sizeof(float), j + 1, in);
      // iterate through the columns of row i
      for (k = 0; k < j + 1; k++){
        // if the correlation value is greater than the given then save the value
        if (k != j && fabs(rowi[k]) > th){
          cutM[params.cutM_index[k] + (used * params.cutM_index[k])] = rowi[k];
        }
      }
    }
    fclose(in);
  }

  // print the cut matrix
  for (i = 0; i < used; i++) {
    for (j = 0; j < used; j++) {
      //printf("%f ", cutM[i * used + j]);
    }
    //printf("\n");
  }

  free(rowi);
  return cutM;
}

/*
 * @param double* l
 * @param int idx1
 * @param int idx2
 */

void swapD(double* l, int idx1, int idx2){
  double temp = l[idx1];
  l[idx1] = l[idx2];
  l[idx2] = temp;
  return;
}

/*
 * @param float* l
 * @param int idx1
 * @param int idx2
 */

void swapF(float* l, int idx1, int idx2){
  float temp = l[idx1];
  l[idx1] = l[idx2];
  l[idx2] = temp;
  return;
}

/*
 * @param double* l
 * @param int size
 */

void quickSortD(double* l, int size){
  if(size<=1) return;
  int pivIdx = (int) size/1.618;//golden ratio
  double pivot = l[pivIdx];
  swapD(l, pivIdx, size-1);
  int leftPlace = 0;
  int i;
  for(i=0;i<size-1;i++){
    if(l[i]<pivot){
      swapD(l, i, leftPlace);
      leftPlace++;
    }
  }
  swapD(l, size-1, leftPlace);
  quickSortD(l,leftPlace);
  quickSortD(&l[leftPlace+1], size-leftPlace-1);
  return;
}

/*
 * @param float* l
 * @param int size
 */

void quickSortF(float* l, int size){
  if(size<=1) return;
  int pivIdx = (int) size/1.618;//golden ratio
  float pivot = l[pivIdx];
  swapF(l, pivIdx, size-1);
  int leftPlace = 0;
  int i;
  for(i=0;i<size-1;i++){
    if(l[i]<pivot){
      swapF(l, i, leftPlace);
      leftPlace++;
    }
  }
  swapF(l, size-1, leftPlace);
  quickSortF(l,leftPlace);
  quickSortF(&l[leftPlace+1], size-leftPlace-1);
  return;
}

/*
 * Calculates the eigenvalues of the given matrix.  This function is a wrapper
 * for the ssyev_ function of the LAPACK package.
 *
 * @param float * smatrix
 *   A pointer to an array of floating point numbers representing the
 *   square n x n similarity matrix.
 * @param int size
 *   The size, n, of the n x n matrix.
 *
 * @return
 *   A pointer to an array of floating point numbers representing the
 *   array of eigenvalues
 */

float* calculateEigen(float * smatrix, int size){

  char jobz = 'N';      // N means don't compute eigenvectors, just eigenvalues
  char uplo = 'U';      // U means the upper matrix is stored
  float * W;            // the array where eignvalues are stored
  float * work;         // a working array. This will be 5 times the size of the final array
  int lwork = 5 * size; // the size of the work array
  int rc;               // indicates the success of the ssyev_ function

  // allocate the arrays
  W    = (float *) malloc(sizeof(float) * size);
  work = (float *) malloc(sizeof(float) * 5 * size);

  ssyev_(&jobz, &uplo , &size, smatrix, &size, W, work, &lwork, &rc);

  // report any errors
  if (rc < 0) {
    printf("\nERROR: During eigenvalue calculation, the %d argument had an illegal value. Continuing anyway...\n", rc);
  }
  else if (rc > 0) {
    printf("\nERROR: The eigenvalue algorithm failed to converge; %d off-diagonal elements of an intermediate tridiagonal form did not converge to zero. Continuing anyway...\n", rc);
  }
  free(work);
  return W;
}

/*
 * returned array will always be sorted and of length size-1
 *
 * @param float* e
 * @param int size
 * @param int m
 */

double * unfolding(float * e, int size, int m){
  int count = 1; // Count equals 1 initially because of 2 lines following loop
                 // which propagates the arrays.
  int i, j = 0;
  double * oX;
  double * oY;
  for(i = 0; i < size - m; i += m) {
    count++;
  }

  oX = (double*) malloc(sizeof(double) * count);
  oY = (double*) malloc(sizeof(double) * count);

  for(i = 0; i < size - m; i += m){
    oX[j] = e[i];
    oY[j] = (i + 1.0) / (double) size;
    j++;
  }
  oX[count-1] = e[size-1];
  oY[count-1] = 1;

  for (i = 1; i < count; i++) {
    if (!(oX[i-1] < oX[i])) {
      printf("\nat postion %d a problem exists\n", i);
      printf("oX[i-1] = %f whilst oX[i] = %f\n", oX[i-1], oX[i]);
    }
  }
  double * yy = (double*) malloc(sizeof(double)*(size));

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, count); //see gsl docs, chapter 27: cspline is a natural spline
  gsl_spline_init(spline, oX, oY, count);

  for (i = 0; i < (size-2); i++) {
    yy[i+1] = gsl_spline_eval(spline, e[i+1], acc);
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  yy[0] = 0.0;
  yy[size-1] = 1.0;
  for (i = 0;i < size - 1; i++) {
    yy[i] = size * (yy[i+1] - yy[i]);
  }
  quickSortD(yy, size-1);
  free(oX);
  free(oY);
  return yy;
}

/**
 * Removes duplicate eigenvalues from an array of eigenvalues.
 *
 * @param float* eigens
 *   The eigenvalue array
 * @param int size
 *   The size of the eigenvalue array
 * @param int* newSize
 *   Upon return, will contain the new size of the updated array
 *
 * @return float *
 *   An array of eigenvalues with duplicates removed.
 */

float* degenerate(float* eigens, int size, int* newSize){
  int i, j = 0;   // used to iterate through the eigenvalues
  int count = 1;  // because one flag is set before the loop
  int * flags;
  float * remDups;

  // iterate through the eigenvalues and change those with a value less
  // that 0.000001 to zero.
  for (i = 0; i < size; i++) {
    if (fabs(eigens[i]) < 0.000001) {
      eigens[i] = 0.0;
    }
  }

  // iterate through the eigenvalues and flag duplicates
  flags = (int*) malloc(sizeof(int) * size);
  memset(flags, 0, size * sizeof(int));
  float temp = eigens[0];
  flags[0] = 1;
  for(i = 1; i < size; i++){
    if(fabs(eigens[i] - temp) > 0.000001){
      count++;
      flags[i] = 1;
      temp = eigens[i];
    }
  }

  // create a new vector without duplicates
  remDups = (float*) malloc(sizeof(float) * count);
  for(i = 0; i < size; i++){
    if(flags[i] == 1){
      remDups[j] = eigens[i];
      j++;
    }
  }
  free(flags);

  // set the newSize argument
  *newSize = count;

  return remDups;
}

/*
 * Returns the averaged Chi-square test across a range of unfolding trials.
 *
 * @param float* eigens
 *   An array of eigenvalues
 * @param int size
 *   The size of the eigenvalue array
 *
 * @return double
 *   A Chi-square value or -1 for failure
 */

double chiSquareTestUnfoldingNNSDWithPoisson(float* eigens, int size, RMTParameters params){
  double chiTest = 0;
  int i = 0;
  int m;

  // We want to generate an average Chi-square value across various levels of
  // unfolding. Therefore, we iterate through the min and max unfolding pace
  // and then average the Chi-square values returned
  for (m = params.minUnfoldingPace; m < params.maxUnfoldingPace; m++) {
    chiTest += chiSquareTestUnfoldingNNSDWithPoisson4(eigens, size, params.nnsdHistogramBin, m, params);

    // if the test failed then return -1
    if (chiTest == -1) {
      return -1;
    }
    i++;
  }

  // return the average Chi-square value
  return chiTest / i;
}

/**
 * Performs a Chi-square test by comparing NNSD of the eigenvalues
 * to
 *
 * @param float* eigens
 *   The eigenvalue array
 * @param int size
 *   The size of the eigenvalue array
 * @param double bin
 *   The relative histogram bin size
 * @param int pace
 *   The unfolding pace
 * @param RMTParameters
 *
 *
 * @return double
 *   A Chi-square value, or -1 on failure
 */

double chiSquareTestUnfoldingNNSDWithPoisson4(float* eigens, int size, double bin, int pace, RMTParameters params){
  int newSize;   // the new size of the eigenvalue array after duplicates removed
  float * newE;  // the new eigenvalue array after duplicates removed
  double * edif;
  double obj;
  double expect;
  double chi = 0;
  int i, j, count;


  // remove duplicates from the list of eigenvalues
  newE = degenerate(eigens, size, &newSize);
  size = newSize;

  // make sure our vector of eigenvalues is still large enough after
  // duplicates have been removed. If not, return a -1
  if (size < params.min_size) {
    printf("    Chi-square test failed: eigenvalue array too small after duplicate removal. See the eigenvector output file.\n");
    return -1;
  }

  edif = unfolding(newE, size, pace);
  free(newE);
  size = size - 1; // see note above unfolding function, will return an array of size-1
  int n = (int) (3.0/bin) + 1;

  for (i = 0; i < n; i++) {
    count = 0;
    for (j=0; j < size; j++) {
      if (edif[j] > i * bin && edif[j] < (i + 1) * bin) {
        count++;
      }
    }
    obj = (double) count;
    expect = (exp(-1 * i * bin) - exp(-1 * (i + 1) *bin)) * size;
    chi += (obj - expect) * (obj - expect) / expect;
  }
  free(edif);
  return chi;
}

/**
 * Prints the command-line usage instructions for the similarity command
 */
void print_threshold_usage() {
  printf("\n");
  printf("Usage: ./RMTGeneNet threshold [options]\n");
  printf("The list of required options:\n");
  printf("  --ematrix|-e The file name that contains the expression matrix.\n");
  printf("                 The rows must be genes or probe sets and columns are samples\n");
  printf("  --method|-m  The correlation method used. Supported methods include\n");
  printf("                 Pearson's correlation and BSpline estimation of Mutual Information.\n");
  printf("                 Provide either 'pc' or mi' as values respectively.\n");
  printf("\n");
  printf("Optional:\n");
  printf("  --th|-t      A decimal indicating the start threshold. For Pearson's.\n");
  printf("                 Correlation (--method pc), the default is 0.99. For Mutual\n");
  printf("                 information (--method mi), the default is the maximum MI value\n");
  printf("                 in the similarity matrix\n");
  printf("  --step|-s    The threshold step size, to subtract at each iteration of RMT.\n");
  printf("                 The default is 0.001\n");
  printf("  --chi|-c     The Chi-square test value which when encountered, RMT will stop.\n");
  printf("                 The default is 200 (corresponds to p-value of 0.01)\n");
  printf("  --perf       Provide this flag to enable performance monitoring.\n");
  printf("\n");
}
