// RandomMatrixModeling.v3.1.c
// uses RMT method to determine proper threshold for correlation matrix
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include "RandomMatrix.h"
#include <time.h>

/**
 * Globals
 */

int* UsedFlag; // holds an array flagging use of indicies in the matrix
int* index1;   // had to rename because index was already taken
int numGenes;
char* inputFileName;
char* inputDir;
int numLinesPerFile;
RMTParameters rmtParameter;
FILE* timing;
// FILE* runInfo;
int verbose; //a boolean for whether or not to print timing data

//set some global params for the RMM stepping
float thresholdStart; //the threshold to start generating cutMatrices with (actually, +thresholdStep)
float thresholdStep;  //the threshold step for each iteration of generation
float chiSoughtValue; //the chiValue the while loop will end on



/**
 * Subroutines
 */

//writes a value to a string, then returns the pointer.  Only supports base 10 numbers.
char* itoa(int val, char* ptr){
  sprintf(ptr,"%d",val);
  return ptr;
} 



int findNumUsed(){
  int i,sum=0;
  for(i=0;i<numGenes;i++){
    if(UsedFlag[i]!=0) sum++;
  }
  return sum;
}

void setIndexArray(){
  int i,j=0;
  for(i=0;i<numGenes;i++){
    if(UsedFlag[i]!=0){
      index1[i]=j;
      j++;
    }
  }
  return;
}





/**
 * The main subroutine
 */
int main(int argc, char** argv) {
  int i, size;

  if(argc==1){
    printf("The arguments for this fucntion are:\n\n");
    printf("'-i': The input file name. Same as used in previous step.\n");
    printf("      Must be the same as the name used in the matrix binary\n");
    printf("      files.  This name will also be used to create files \n");
    printf("      associated with this run of the program.\n\n");
    printf("Optional:\n\n");
    printf("'-b': The initial threshold(+1*step) value that will be used.\n");
    printf("      [b=Begin value] Default: 0.9200\n");
    printf("'-s': The threshold step size used each iteration. Default: 0.001\n");
    printf("'-c': The chi-square test value that the loop will stop on.\n");
    printf("      Default: 200\n\n");
    printf("'-v': Set the performance collection. Has two values possible values,\n");
    printf("      ON/OFF . [v=Verbose] Default: ON\n");
    printf("Examples:\n");
    printf("<executable> -i <input.file.name> \n");
    printf("<exec> -i <input.file.name> -s <0.0001> -v ON\n\n");
    printf("Note: Order is not important, but spaces are required between the \n");
    printf("flag and its value.  Also note that erroneous inputs could cause the\n");
    printf("program to crash without a helpful error message...so get it right.\n\n");
    printf("Exiting.\n");
    exit(0);
  }

  inputDir = "Pearson/";

  // initialize default values, which may be overwritten in the command
  verbose=1;
  thresholdStart = 0.92;
  thresholdStep = 0.001;
  chiSoughtValue = 200;

  // parse the command line flags
  for (i = 1; i < argc; i += 2) {
    if (argv[i][1] == 'i'){
      inputFileName = argv[i+1];
    }
    else if (argv[i][1] == 'b') {
      thresholdStart = atof(argv[i+1]);
    }
    else if (argv[i][1] == 's') {
      thresholdStep = atof(argv[i+1]);
    }
    else if (argv[i][1] == 'c') {
      chiSoughtValue = atof(argv[i+1]);
    }
    else if (argv[i][1] == 'v') {
      if (strcmp(argv[i+1], "ON") == 0 || strcmp(argv[i+1], "on") == 0 || strcmp(argv[i+1],"On") == 0 || strcmp(argv[i+1],"oN") == 0) {
        verbose = 1;
      }
      else{
        verbose = 0;
      }
    }
    else {
      printf("Flag '%s' not recognized, exiting.", argv[i]);
      exit(0);
    }
  }

  // read from first matrix file to determine some stuff
  char* filename;
  char num[4];
  int len = strlen(inputFileName);
  len += strlen(inputDir);
  len += strlen(".pc");
  len += strlen("xxx.bin");
  len++;
  filename = (char*)malloc(sizeof(char)*len);
  memset(filename, '\0', len);
  strcat(filename, inputDir);
  strcat(filename, inputFileName);
  strcat(filename, ".pc");
  strcat(filename, itoa(0, num));
  strcat(filename, ".bin");
  FILE* info;
  info = fopen(filename, "rb");
  fread(&numGenes, sizeof(int), 1, info);
  fread(&numLinesPerFile, sizeof(int), 1, info);
  fclose(info);
  free(filename);

  // initialize some global file pointers that will be used throughout
  // runInfo = fileOpenHelper(".runInfo.txt");

  if(verbose == 1){
    timing = fileOpenHelper(".timing.txt");
  }

  // print some preliminary run info to the runInfo file
  // fprintf(runInfo, "number of genes: %d\n", numGenes);
  // fprintf(runInfo, "threshold start: %f\n", thresholdStart);
  // fprintf(runInfo, "threshold step:  %f\n", thresholdStep);
  // fprintf(runInfo, "sought chi val:  %f\n\n", chiSoughtValue);
  // fprintf(runInfo, "Below are the thresholds and corresponding cut matrix sizes\n");

  //initialize the global RMTParameters struct
  rmtParameter.nnsdHistogramBin = 0.05;
  rmtParameter.chiSquareTestThreshold = 99.607;
  rmtParameter.minUnfoldingPace = 10;
  rmtParameter.maxUnfoldingPace = 41;
  rmtParameter.mimiModuleSize = 4;
  rmtParameter.edHistogramBin = 0.1;

  //allocate memory for UsedFlags
  UsedFlag = (int *) malloc(numGenes * sizeof(int));

  //allocate memory for index
  index1 = (int *) malloc(numGenes * sizeof(int));

  time_t start, end;
  fflush(timing);
  time(&start);
  int rc = determinePearsonCorrelationThreshold_LargeMatrix();
  time(&end);
  if (verbose) {
    fprintf(timing, "Minutes in determining: %f\n", (end - start) / 60.0);
  }
  free(index1);
  free(UsedFlag);

  // fclose(runInfo);
  if(verbose){
    fclose(timing);
  }
  exit(rc);
}


/**
 *
 */
int determinePearsonCorrelationThreshold_LargeMatrix() {
  float* newM;
  int size;
  time_t start, end;
  FILE* eigenF, *chiF;

  // open the output files
  // eigenF = fileOpenHelper(".eigens.txt");
  chiF = fileOpenHelper(".chiVals.txt");
  fprintf(chiF, "Threshold\tChi-square\tCut Matrix Size\n");  

  float th = thresholdStart;
  double finalTH  = 0.0;
  double finalChi = 10000.0;
  double minTH    = 1.0;
  double minChi   = 10000.0;
  double maxTH    = 0.0;
  double maxChi   = 0.0;
  int i           = 0;
  double chi;
  float* E;  //array for eigenvalues

  do {
    // decrement the threshold using the step value and then retreive the 
    // matrix that contains only the threshold and above
    th = th - thresholdStep;
    printf("Testing threshold: %f...\n", th);
//    time(&start);
    printf("  reading bin files...\n");
    newM = readPearsonCorrelationMatrix(th, &size);
//    time(&end);
//    if(verbose){
//      fprintf(timing,"Minutes for cut matrix: %f\n", (end-start)/60.0);
//    }
//    fprintf(runInfo,"%f\t%d\n",th,size);
//    fflush(runInfo);
    if (size >= 100) {
//      time(&start);
      printf("  calculating eigenvalues for n x n matrix of size n = %d...\n", size);
      E = calculateEigen(newM, size);
//      time(&end);

      free(newM);
//      if(verbose){
//        fprintf(timing,"Minutes calculating Eigens: %f\n", (end-start)/60.0);
//      }
      /* print out eigenvalues to file */
      // fprintf(eigenF, "%f\t", th);
      // for(i=0 ; i<size ; i++){
        // fprintf(eigenF, "%f\t", E[i]);
      // }
      // fprintf(eigenF,"\n");
//      time(&start);
      printf("  testing similarity of NNSD with Poisson...\n");
      chi = chiSquareTestUnfoldingNNSDWithPoisson(E, size, rmtParameter.nnsdHistogramBin, rmtParameter.minUnfoldingPace, rmtParameter.maxUnfoldingPace);
//      time(&end);
//      if(verbose){
//        fprintf(timing,"Minutes Chi Square Testing: %f\n", (end-start)/60.0);
//      }
//      fflush(timing);
      fprintf(chiF, "%f\t%f\t%d\n", th, chi, size);
      fflush(chiF);
      free(E);
      printf("  chi = %f\n", chi);
      if(chi < minChi){
        minChi = chi;
        minTH = th;
      }
      if (chi < rmtParameter.chiSquareTestThreshold){
        finalTH = th;
        finalChi = chi;
      }
      if (finalChi < rmtParameter.chiSquareTestThreshold && chi > finalChi && th < finalTH){
        maxChi = chi;
        maxTH = th;
      }
    }
    else{
      free(newM);
    }
  }
  while(maxChi < chiSoughtValue || size == numGenes);

  /*
   If finalChi is still greater than threshold, check the small scale
  */

  if (finalChi > rmtParameter.chiSquareTestThreshold){
  // fprintf(runInfo,"checking small scale\n");
    fprintf(chiF, "checking small scale\n");
    th = (float)minTH + 0.2;
    for(i = 0 ; i <= 40 ; i++){
      th = th - thresholdStep*i;
//      time(&start);
      newM = readPearsonCorrelationMatrix(th, &size);
//      time(&end);
//      if(verbose){
//        fprintf(timing,"Minutes for cut matrix: %f\n", (end-start)/60.0);
//      }
      // fprintf(runInfo,"%f\t%d\n",th,size);
      if(size >= 100){
//        time(&start);
        E = calculateEigen(newM, size);
//        time(&end);
//        if(verbose){
//          fprintf(timing,"Minutes for cut matrix: %f\n", (end-start)/60.0);
//        }
        free(newM);
        /* print out eigenvalues to file */
        // fprintf(eigenF, "%f\t", th);
        // for(i=0 ; i<size ; i++){
          // fprintf(eigenF, "%f\t", E[i]);
        // }
        // fprintf(eigenF,"\n");
//        time(&start);
        chi = chiSquareTestUnfoldingNNSDWithPoisson(E, size, rmtParameter.nnsdHistogramBin, rmtParameter.minUnfoldingPace, rmtParameter.maxUnfoldingPace);
//        time(&end);
//        if(verbose){
//          fprintf(timing,"Minutes Chi Square Testing: %f\n", (end-start)/60.0);
//        }
        fprintf(chiF, "%f\t%f\t%d\n", th, chi, size);
        free(E);
        if(chi < minChi){
          minChi = chi;
          minTH = th;
        }
        if (chi < rmtParameter.chiSquareTestThreshold){
          finalTH = th;
          finalChi = chi;
        }
      }//end if size >= 100
      else{
        free(newM);
      }
    }//end for 1 -> 40 loop
  }//end if finalChi > rmt...

  //close the chi and eigen files now that results are written
  fclose(chiF);
  // fclose(eigenF);
  /*
  Set the Properties file according to success or failure
  */
  // fprintf(runInfo, "==================================\n");
  if(finalChi < rmtParameter.chiSquareTestThreshold){
    finalTH = ceil(finalTH * 10000) / 10000.0;
    FILE* th;
    th = fileOpenHelper(".th.txt");
    fprintf(th, "%f", finalTH);
    fclose(th);
    return 0;
  }
  else{
    finalTH = ceil(finalTH * 10000) / 10000.0;
    return -2;
  }
}

/**
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
 *
 */
float * readPearsonCorrelationMatrix(float th, int * size) {

  float * rowj;
  char * filename; // the name of the correlation matrix bin file
  char num[4];
  int len;         // the length of the filename and path to the correlation bin file
  int i, h;        // used to iterate through the bin files
  int j;           // used to iterate through the rows of each bin file
  int k;           // used to iterate through the cols of each bine file
  int z;           // the number of binary files
  int used;        // holds the number of genes (probesets) that have a greater thrshold
  int limit;       // the maximum row in the current correlation bin file
  int junk;        // a dummy variable
  FILE* in;

  memset(index1, -1, sizeof(int) * (numGenes));

  rowj = (float*) malloc(sizeof(float) * numGenes);

  // reserve the proper amount of memory for the filename and path
  len =  strlen(inputFileName);
  len += strlen(inputDir);
  len += strlen(".pc");
  len += strlen("xxx.bin");
  len++;
  filename = (char*) malloc(sizeof(char) * len);
  memset(filename, '\0', len);

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
  z = (numGenes - 1) / numLinesPerFile;
  for (i = 0; i <= z; i++) {
    memset(filename, '\0', len);
    strcat(filename, inputDir);
    strcat(filename, inputFileName);
    strcat(filename, ".pc");
    strcat(filename, itoa(i, num));
    strcat(filename, ".bin");
    in = fopen(filename, "rb");
    fread(&junk, sizeof(int), 1, in); // numGenes
    fread(&junk, sizeof(int), 1, in); // numLinesPerFile
    if (i != z) {
      limit = (i + 1) * numLinesPerFile;
    }
    else{
      limit = numGenes;
    }

    // iterate through the rows and columns of the file and look for 
    // entries greater than the provided threshold.  When found, use
    // the row and column indexes to set a '1' in the UsedFlag array.
    // this array indicates which genes have values we want to keep. 
    for (j = i * numLinesPerFile; j < limit; j++) {
      fread(rowj, sizeof(float), j + 1, in);
      for (k = 0; k < j + 1; k++) {
        // if the correlation value is greater than the given threshold then 
        // flag the row/column indexes
        if (k != j && fabs(rowj[k]) > th) {
          UsedFlag[k] = 1;
          UsedFlag[j] = 1;
        }
      }
    }
    fclose(in);
  }

  // get the number of genes (or probesets) that have a correlation value
  // greater than the provided threshold value
  used = findNumUsed();
  setIndexArray();

  // now that we know how many genes have a threshold greater than the
  // given we can allocate memory for new cut matrix
  float* cutM = (float*)calloc(used*used, sizeof(float));

  // initialize the diagonal to 1
  for (i = 0; i < used; i++) {
    cutM[i + i * used] = 1;
  }

  // set the incoming size argument to be the size dimension of the cut matrix
  *size = used;

  // ------------------------------
  // Step #2: Now build the cut matrix by retreiving the correlation values
  // for each of the genes identified previously
  for(h = 0; h < z; h++){
    memset(filename,'\0', len );
    strcat(filename, inputDir);
    strcat(filename, inputFileName);
    strcat(filename, ".pc");
    strcat(filename, itoa(h, num));
    strcat(filename, ".bin");
    in = fopen(filename, "rb");
    fread(&junk, sizeof(int), 1, in); // numGenes
    fread(&junk, sizeof(int), 1, in); // numLinesPerFile
    if (i != z) {
      limit = (h + 1) * numLinesPerFile;
    }
    else{
      limit = numGenes;
    }

    // iterate through the rows of the bin file
    for (i = h * numLinesPerFile; i < limit; i++) {
      fread(rowj, sizeof(float), i+1, in);
      // iterate through the columns of the bin file
      for (j = 0; j < i + 1; j++){
        // if the correlation value is greater than the given then save the value
        if (i != j && fabs(rowj[j]) > th){
          cutM[index1[j] + used * index1[i]] = rowj[j];
        }
      }
    }
    fclose(in);
  }

  free(rowj);
  free(filename);
  return cutM;
}

/**
 * Opens a file using the input file name as theprefix and a given string
 * as a suffix.
 *
 * @param char * extension
 *   A string to be used as the file suffix.
 *
 * @return 
 *   A file pointer
 */
FILE * fileOpenHelper(char* extension) {
  char * filename; // filename buffer
  FILE * fp;       // file pointer to be returned
  int len;         // holds the length of the newly created file name

  // determine the size of the file name and build the filename
  len = strlen(inputFileName);
  len += strlen(extension);
  len++;
  filename = (char *) malloc(sizeof(char) * len);
  memset(filename,'\0', len);
  strcat(filename,inputFileName);
  strcat(filename,extension);

  // open the file
  fp = fopen(filename, "w");
  if (fp == NULL) {
    printf("\nUnable to create file '%s' . Please check input filename for validity.\nExiting.\n", filename);
    exit(0);
  }
  free(filename);
  return fp;
}
