#include "RMTThreshold.h"

RMTThreshold::RMTThreshold(EMatrix * ematrix, char ** method, int num_methods,
    char * th_method, float thresholdStart, float thresholdStep, float chiSoughtValue,
    char * clustering, int min_cluster_size, int max_missing, int max_modes,
    float min_range)
  : ThresholdMethod(ematrix, method, num_methods, th_method, clustering,
      min_cluster_size, max_missing, max_modes, min_range) {

  this->thresholdStart = thresholdStart;
  this->thresholdStep  = thresholdStep;
  this->chiSoughtValue = chiSoughtValue;
  this->clustering = clustering;

  //
  minEigenVectorSize = 50;

  // TODO: perhaps the user should have a bit more control over what these
  // values should be? 99.607 (Chi-square df = 60 (number of bins in NNSD
  // histogram), p-value = 0.001). What if user wants a p-value of 0.01?
  nnsdHistogramBin       = 0.05;
  chiSquareTestThreshold = 99.607;
  minUnfoldingPace       = 10;
  maxUnfoldingPace       = 41;

  finalTH  = 0.0;
  finalChi = 10000.0;
  minTH    = 1.0;
  minChi   = 10000.0;
  maxChi   = 0.0;
}

/**
 *
 */
RMTThreshold::~RMTThreshold() {

}
/*
 *
 *
 */
float RMTThreshold::findThreshold() {

  // A reduced expression matrix containing only gene-paris with
  // correlations >= the test threshold.
  float* newM;
  // The size of newM.
  int size;
  // File handles for the eigenvector and Chi-square output files.
  FILE *eigenF;
  FILE *chiF;
  // The output file name
  char chi_filename[1024];
  // The output file name
  char eigen_filename[1024];
  // The threshold currently being tested.
  float th = thresholdStart;
  // The current chi-square value for the threshold being tested.
  float chi;
  // Array for eigenvalues.
  float * E;
  // The output file prefix.
  char * file_prefix = ematrix->getFilePrefix();
  // The number of samples in the expression matrix.
  int num_samples = ematrix->getNumSamples();


  // Open the output files and print the headers.
  if (clustering) {
    if (max_missing > num_samples) {
      sprintf(chi_filename, "%s.%s.mcs%d.md%d.mmINF.chiVals.txt", file_prefix, th_method, min_cluster_size, max_modes);
      sprintf(eigen_filename, "%s.%s.mcs%d.md%d.mmINF.eigenVals.txt", file_prefix, th_method, min_cluster_size, max_modes);
    }
    else {
      sprintf(chi_filename, "%s.%s.mcs%d.md%d.mm%d.chiVals.txt", file_prefix, th_method, min_cluster_size, max_modes, max_missing);
      sprintf(eigen_filename, "%s.%s.mcs%d.md%d.mm%d.eigenVals.txt", file_prefix, th_method, min_cluster_size, max_modes, max_missing);
    }
  }
  else {
    sprintf(chi_filename, "%s.%s.chiVals.txt", file_prefix, th_method);
    sprintf(eigen_filename, "%s.%s.eigenVals.txt", file_prefix, th_method);
  }


  chiF = fopen(chi_filename, "w");
  eigenF = fopen(eigen_filename, "w");
  fprintf(chiF, "Threshold\tChi-square\tCut Matrix Size\n");

  // Iterate through successively smaller threshold values until the following
  // conditions are met:
  // 1)  A Chi-square of
  do {
    printf("\n");
    printf("  testing threshold: %f...\n", th);

    if (!clustering) {
      newM = read_similarity_matrix_bin_file(th, &size);
    }
    else {
      newM = read_similarity_matrix_cluster_file(th, &size);
    }

    printf("  found matrix of size n x n, n = %d...\n", size);
    if (size >= minEigenVectorSize) {
      printf("  calculating eigenvalues...\n");
      E = calculateEigen(newM, size);
      free(newM);

      // print out eigenvalues to file
      fprintf(eigenF, "%f\t", th);
      for (int i = 0; i < size ; i++) {
        fprintf(eigenF, "%f\t", E[i]);
      }
      fprintf(eigenF,"\n");
      fflush(eigenF);

      printf("  testing similarity of NNSD with Poisson...\n");
      chi = getNNSDChiSquare(E, size);
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
        if (chi < chiSquareTestThreshold){
          finalTH = th;
          finalChi = chi;
        }
        if (finalChi < chiSquareTestThreshold && chi > finalChi && th < finalTH){
          maxChi = chi;
        }
      }
    }
    else{
      free(newM);
    }

    // decrement the threshold using the step value and then retrieve the
    // matrix that contains only the threshold and higher
    th = th - thresholdStep;
  }
  while(maxChi < chiSoughtValue);


  // If finalChi is still greater than threshold, check the small scale
  if (finalChi > chiSquareTestThreshold) {
    fprintf(chiF, "checking small scale\n");
    th = (float)minTH + 0.2;
    for (int i = 0 ; i <= 40 ; i++) {
      th = th - thresholdStep * i;
      newM = read_similarity_matrix_cluster_file(th, &size);

      if (size >= 100) {
        E = calculateEigen(newM, size);
        free(newM);

        // print out eigenvalues to file
        fprintf(eigenF, "%f\t", th);
        for (int j = 0 ; j < size ; j++) {
          fprintf(eigenF, "%f\t", E[j]);
        }
        fprintf(eigenF, "\n");
        chi = getNNSDChiSquare(E, size);
        fprintf(chiF, "%f\t%f\t%d\n", th, chi, size);
        fflush(chiF);
        free(E);

        if (chi < minChi) {
          minChi = chi;
          minTH = th;
        }
        if (chi < chiSquareTestThreshold) {
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
  if(finalChi < chiSquareTestThreshold){
    finalTH = ceil(finalTH * 10000) / 10000.0;
    FILE* th;
    char filename[1024];
    if (clustering) {
      if (max_missing > num_samples) {
        sprintf(filename, "%s.%s.mcs%d.md%d.mmINF.th.txt", file_prefix, th_method, min_cluster_size, max_modes);
      }
      else {
        sprintf(filename, "%s.%s.mcs%d.md%d.mm%d.th.txt", file_prefix, th_method, min_cluster_size, max_modes, max_missing);
      }
    }
    else {
      sprintf(filename, "%s.%s.th.txt", file_prefix, th_method);
    }
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
 * Parses the similarity matrix stored in the binary file format.
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

float * RMTThreshold::read_similarity_matrix_bin_file(float th, int * size) {

  float * cutM;    // the resulting cut similarity matrix
  float * rowj;    // holds the float value from row i in the bin file
  int i;           // used to iterate through the bin files
  int j;           // used to iterate through the rows of each bin file
  int k;           // used to iterate through the cols of each bine file
  int z;           // the number of binary files
  int used;        // holds the number of genes (probesets) that have a greater thrshold
  int limit;       // the maximum row in the current correlation bin file
  int junk;        // a dummy variable
  FILE* in;
  char filename[1024]; // used for storing the bin file name

  int file_num_genes;
  int file_num_lines;


  // open the file and get the number of genes and the lines per file
  // these data are the first two integers in the file
  FILE* info;
  sprintf(filename, "%s/%s.%s%d.bin", bin_dir, ematrix->getFilePrefix(), th_method, 0);
  // TODO: check that file exists before trying to open
  info = fopen(filename, "rb");
  fread(&file_num_genes, sizeof(int), 1, info);
  fread(&file_num_lines, sizeof(int), 1, info);
  fclose(info);

  printf("  Genes: %d, Num lines per file: %d\n", file_num_genes, file_num_lines);

  // Initialize vector for holding which gene clusters will be used. This will
  // tell us how big the cut matrix needs to be.
  int num_genes = ematrix->getNumGenes();
  int * usedFlag = (int *) malloc(num_genes * sizeof(int));
  int * cutM_index = (int *) malloc(num_genes * sizeof(int));

  memset(cutM_index, -1, sizeof(int) * (num_genes));
  memset(usedFlag, 0, sizeof(int) * (num_genes));


  rowj = (float *) malloc(sizeof(float) * file_num_genes);

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
  z = (file_num_genes - 1) / file_num_lines;
  for (i = 0; i <= z; i++) {

    sprintf(filename, "%s/%s.%s%d.bin", bin_dir, ematrix->getFilePrefix(), th_method, i);
    in = fopen(filename, "rb");
    fread(&junk, sizeof(int), 1, in); // file_num_genes
    fread(&junk, sizeof(int), 1, in); // file_num_lines
    if (i != z) {
      limit = (i + 1) * file_num_lines;
    }
    else{
      limit = file_num_genes;
    }

    // iterate through the rows and columns of the file and look for
    // entries greater than the provided threshold.  When found, use
    // the row and column indexes to set a '1' in the  array.
    // this array indicates which genes have values we want to keep.
    for (j = i * file_num_lines; j < limit; j++) {
      fread(rowj, sizeof(float), j + 1, in);
      for (k = 0; k < j + 1; k++) {
        // if the correlation value is greater than the given threshold then
        // flag the row/column indexes
        if (k != j && fabs(rowj[k]) > th) {
          usedFlag[k] = 1;
          usedFlag[j] = 1;
        }
      }
    }
    fclose(in);
  }

  // get the number of genes (or probe sets) that have a correlation value
  // greater than the provided threshold value
  used = 0;
  j = 0;
  for (i = 0; i < file_num_genes; i++) {
    if (usedFlag[i] == 1) {
      used++;
      cutM_index[i] = j;
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

    sprintf(filename, "%s/%s.%s%d.bin", bin_dir, ematrix->getFilePrefix(), th_method, i);
    in = fopen(filename, "rb");
    fread(&junk, sizeof(int), 1, in); // file_num_genes
    fread(&junk, sizeof(int), 1, in); // file_num_lines
    if (i != z) {
      limit = (i + 1) * file_num_lines;
    }
    else{
      limit = file_num_genes;
    }
    // iterate through the rows of the bin file
    for (j = i * file_num_lines; j < limit; j++) {
      fread(rowj, sizeof(float), j + 1, in);
      // iterate through the columns of row i
      for (k = 0; k < j + 1; k++){
        // if the correlation value is greater than the given then save the value
        if (k != j && fabs(rowj[k]) > th){
          cutM[cutM_index[k] + (used * cutM_index[j])] = rowj[k];
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

  // free memory
  free(rowj);
  free(cutM_index);
  free(usedFlag);
  return cutM;
}

/*
 * Parses the similarity matrix stored in the clusters file format.
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

float * RMTThreshold::read_similarity_matrix_cluster_file(float th, int * size) {

  // The resulting cut similarity matrix.
  float * cutM;

  // The maximum number of clusters that any gene pair is allowed to have.
  int max_clusters = 15;

  // Initialize vector for holding which gene clusters will be used. This will
  // tell us how big the cut matrix needs to be.
  int num_genes = ematrix->getNumGenes();
  int num_samples = ematrix->getNumSamples();
  int * usedFlag = (int *) malloc(num_genes * max_clusters * sizeof(int));
  int * cutM_index = (int *) malloc(num_genes * max_clusters * sizeof(int));
  const char * filename;

  memset(cutM_index, -1, sizeof(int) * (num_genes * max_clusters));
  memset(usedFlag, 0, sizeof(int) * (num_genes * max_clusters));

  // Make sure the output directory exists.
  struct stat st = {0};
  char clusterdir[100];
  sprintf(clusterdir, "clusters-%s", th_method);
  if (stat(clusterdir, &st) == -1) {
    fprintf(stderr, "The clusters directory is missing. Cannot continue.\n");
    exit(-1);
  }

  // ------------------------------
  // Step #1: iterate through the cluster files to find out how many
  // genes there will be in the cut matrix.
  int limit = (int) (th * 100);
  int i;
  for (i = 100; i >= limit; i--) {
    char dirname[1024];
    sprintf(dirname, "./%s/%03d", clusterdir, i);

    DIR * dir;
    dir = opendir(dirname);
    if (dir) {
      struct dirent * entry;
      while ((entry = readdir(dir)) != NULL) {
        filename = entry->d_name;

        // Skip the . and .. files.
        if (strcmp(filename, ".") == 0 || strcmp(filename, "..") == 0) {
          continue;
        }

        // The file must end in a suffix with 3 digits followed by .txt.
        // We use a regular expression to see if this is true. If not, then
        // skip this file.
        regex_t re;
        char pattern[64] = "[0-9][0-9][0-9].txt";
        int  rc;
        size_t nmatch = 2;
        regmatch_t pmatch[2];
        if (0 != (rc = regcomp(&re, pattern, 0))) {
          printf("regcomp() failed, returning nonzero (%d)\n", rc);
          exit(-1);
        }
        if (regexec(&re, filename, nmatch, pmatch, 0)) {
          regfree(&re);
          continue;
        }
        regfree(&re);

        // Construct the full path to the file.
        char path[1024];
        sprintf(path, "%s/%s", dirname, filename);

        // Open each file and traverse the elements
        FILE * fp = fopen(path, "r");
        if (!fp) {
          fprintf(stderr, "Can't open file, %s. Cannot continue.\n", path);
          exit(-1);
        }
        int j, k, cluster_num, num_clusters, cluster_samples, num_missing, num_outliers, num_goutliers, num_threshold;
        char samples[num_samples];
        char cscores[255];
        float cv;
        while (!feof(fp)) {

          // Read in the fields for this line. We must read in 11 fields or
          // we will skip the line.
          int matches = fscanf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s", &j, &k, &cluster_num, &num_clusters, &cluster_samples, &num_missing, &num_outliers, &num_goutliers, &num_threshold, (char *) &cscores, (char *) &samples);
          if (matches < 11) {
            char tmp[num_samples*2];
            matches = fscanf(fp, "%s\n", (char *)&tmp);
            continue;
          }

          // Get the score for the selected method.
          float ** scores = parseScores((char *) &cscores);
          cv = *(scores[th_method_index]);

          // Make sure there aren't more clusters than was allowed.
          if (cluster_num > max_clusters) {
            fprintf(stderr, "Currently, only %d clusters are supported. Gene pair (%i, %i) has %d clusters.\n", max_clusters, j, k, cluster_num);
            exit(-1);
          }

          // Get the range of expression for this cluster.
          float range_j = this->ematrix->getRange(j-1, samples);
          float range_k = this->ematrix->getRange(k-1, samples);
          int range_okay = 1;
          if (min_range > 0 && (range_j <= min_range || range_k <= min_range)) {
            range_okay = 0;
          }

          // filter the record.
          if (fabs(cv) >= th && cluster_samples >= min_cluster_size  &&
              num_missing <= max_missing && num_clusters <= max_modes &&
              range_okay == 1) {

            usedFlag[(j - 1) * max_clusters + (cluster_num - 1)] = 1;
            usedFlag[(k - 1) * max_clusters + (cluster_num - 1)] = 1;
          }

          for (int l = 0; l < num_methods; l++) {
            free(scores[l]);
          }
          free(scores);
        }
        fclose(fp);
      }
      closedir(dir);
    }
    else {
      fprintf(stderr, "The clusters sub directory, %s, is missing. Cannot continue.\n", dirname);
      exit(-1);
    }
  }

  // Count the number of used genes and store in the cutM_index array
  // where each gene cluster will go in the cut matrix.
  int used = 0;
  int j = 0;
  for (i = 0; i < num_genes * max_clusters; i++) {
    if (usedFlag[i] == 1) {
      used++;
      cutM_index[i] = j;
      j++;
    }
  }

  // now that we know how many genes have a threshold greater than the
  // given we can allocate memory for new cut matrix
  cutM = (float *) malloc(sizeof(float) * used * used);
  memset(cutM, 0, sizeof(float) * used * used);
  // initialize the diagonal to 1
  for (i = 0; i < used; i++) {
    cutM[i + i * used] = 1;
  }

  // set the incoming size argument to be the size dimension of the cut matrix
  *size = used;

  // ------------------------------
  // Step #2: Now build the cut matrix by retrieving the correlation values
  // for each of the genes identified previously.
  for (i = 100; i >= limit; i--) {
    char dirname[1024];
    sprintf(dirname, "./%s/%03d", clusterdir, i);

    DIR * dir;
    dir = opendir(dirname);
    if (!dir) {
      fprintf(stderr, "The clusters sub directory, %s, is missing. Cannot continue.\n", dirname);
      continue;
    }
    struct dirent * entry;
    while ((entry = readdir(dir)) != NULL) {
      const char * filename = entry->d_name;

      // Skip the . and .. files.
      if (strcmp(filename, ".") == 0 || strcmp(filename, "..") == 0) {
        continue;
      }

      // The file must end in a suffix with 3 digits followed by .txt.
      // We use a regular expression to see if this is true. If not, then
      // skip this file.
      regex_t re;
      char pattern[64] = "[0-9][0-9][0-9].txt";
      int  rc;
      size_t nmatch = 2;
      regmatch_t pmatch[2];
      if (0 != (rc = regcomp(&re, pattern, 0))) {
        printf("regcomp() failed, returning nonzero (%d)\n", rc);
        exit(-1);
      }
      if (regexec(&re, filename, nmatch, pmatch, 0)) {
        regfree(&re);
        continue;
      }
      regfree(&re);

      // Construct the full path to the file.
      char path[1024];
      sprintf(path, "%s/%s", dirname, filename);

      // Open each file and traverse the elements
      FILE * fp = fopen(path, "r");
      if (!fp) {
        fprintf(stderr, "Can't open file, %s. Cannot continue.\n", path);
        exit(-1);
      }
      int j, k, cluster_num, num_clusters, cluster_samples, num_missing, num_outliers, num_goutliers, num_threshold;
      char samples[num_samples];
      char cscores[255];
      float cv;
      while (!feof(fp)) {

        // Read in the fields for this line. We must read in 11 fields or
        // we will skip the line.
        int matches = fscanf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s", &j, &k, &cluster_num, &num_clusters, &cluster_samples, &num_missing, &num_outliers, &num_goutliers, &num_threshold, (char *) &cscores, (char *) &samples);
        if (matches < 11) {
          char tmp[num_samples*2];
          matches = fscanf(fp, "%s\n", (char *)&tmp);
          continue;
        }

        // Get the score for the selected method.
        float ** scores = parseScores((char *) &cscores);
        cv = *(scores[th_method_index]);

        if (cluster_num > max_clusters) {
          fprintf(stderr, "Currently, only %d clusters are supported. Gene pair (%i, %i) has %d clusters.\n", max_clusters, j, k, cluster_num);
          exit(-1);
        }

        // Get the range of expression for this cluster.
        float range_j = this->ematrix->getRange(j-1, samples);
        float range_k = this->ematrix->getRange(k-1, samples);
        int range_okay = 1;
        if (min_range > 0 && (range_j <= min_range || range_k <= min_range)) {
          range_okay = 0;
        }

        if (fabs(cv) >= th && cluster_samples >= min_cluster_size  &&
            num_missing <= max_missing && num_clusters <= max_modes &&
            range_okay == 1) {

          int index_j = cutM_index[(j - 1) * max_clusters + (cluster_num - 1)];
          int index_k = cutM_index[(k - 1) * max_clusters + (cluster_num - 1)];
          int ci = index_k + (used * index_j);
          cutM[ci] = cv;
        }
      }
      fclose(fp);
    }
    closedir(dir);
  }

  // free memory
  free(cutM_index);
  free(usedFlag);
  return cutM;
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

float * RMTThreshold::calculateEigen(float * smatrix, int size) {

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

/**
 *
 * @param float* e
 *   A sorted eigenvalue array with duplicates removed.
 * @param int size
 *   The size of the eigenvalue array.
 * @param int m
 *  The pace.
 *
 * @return
 *  Returned array will always be sorted and of length size-1
 */

float * RMTThreshold::unfolding(float * e, int size, int m) {
  // Count equals 1 initially because of 2 lines following loop
  // which propagates the arrays.
  int count = 1;
  int i, j = 0;

  // Figure out how many points we will use from the submitted
  // eigenvalue array.  If the pace (m) is 10 and the size is 100
  // then the count will be 10 and we will use 10 points for spline
  // calculation
  for(i = 0; i < size - m; i += m) {
    count++;
  }

  // Retrieve the 'count' number of points from the eigenvalue array.
  double *oX = (double*) malloc(sizeof(double) * count);
  double *oY = (double*) malloc(sizeof(double) * count);
  for(i = 0; i < size - m; i += m){
    oX[j] = e[i];
    oY[j] = (i + 1.0) / (double) size;
    j++;
  }
  oX[count-1] = e[size-1];
  oY[count-1] = 1;

  // Make sure all points are in increasing order. If not then a problem.
  for (i = 1; i < count; i++) {
    if (!(oX[i-1] < oX[i])) {
      printf("\nat postion %d a problem exists\n", i);
      printf("oX[i-1] = %f whilst oX[i] = %f\n", oX[i-1], oX[i]);
    }
  }

  // Initialize the spline function using a csspline. See gsl docs,
  // chapter 27: cspline is a natural spline.
  // Changed to akima spline (2/25/2016) -- spf
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
//  gsl_spline *spline = gsl_spline_alloc(gsl_interp_akima, count);
//  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, count);
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_akima_periodic, count);
  gsl_spline_init(spline, oX, oY, count);

  // Estimate new eigenvalues along the spline curve.  We provide the
  // eigenvalues and the interplolation function gives us new y values
  // in the range of 0 to 1.
  float * yy = (float*) malloc(sizeof(float) * size);
  for (i = 1; i < size - 1; i++) {
    yy[i] = (float) gsl_spline_eval(spline, e[i], acc);
  }

  // Calculate the spacing array.
  yy[0] = 0.0;
  yy[size - 1] = 1.0;
  for (i = 0; i < size - 1; i++) {
    yy[i] = (yy[i+1] - yy[i]) * size;
  }
  quickSortF(yy, size-1);

  // Free up memory.
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  free(oX);
  free(oY);

  // Return the nearest neighbor spacing array.
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
 */#include <getopt.h>

float * RMTThreshold::degenerate(float* eigens, int size, int* newSize) {
  int i, j = 0;   // used to iterate through the eigenvalues
  int count = 1;  // because one flag is set before the loop
  int * flags;
  float * remDups;

  // Iterate through the eigenvalues and change those with a value less
  // that 0.000001 to zero.
  for (i = 0; i < size; i++) {
    if (fabs(eigens[i]) < 0.000001) {
      eigens[i] = 0.0;
    }
  }

  // Iterate through the eigenvalues and flag duplicates.
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

  // Create a new vector without duplicates.
  remDups = (float*) malloc(sizeof(float) * count);
  for(i = 0; i < size; i++){
    if(flags[i] == 1){
      remDups[j] = eigens[i];
      j++;
    }
  }
  free(flags);


  // Set the newSize argument.
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
 * @return float
 *   A Chi-square value or -1 for failure
 */

float RMTThreshold::getNNSDChiSquare(float* eigens, int size) {
  float chiTest = 0;
  float avg_chiTest;
  int i = 0;
  int m;
  // The new eigenvalue array after duplicates removed
  float * newE;
  // The new size of the eigenvalue array after duplicates removed
  int newSize;

  // Remove duplicates from the list of eigenvalues.
  newE = degenerate(eigens, size, &newSize);

  // Make sure our vector of eigenvalues is still large enough after
  // duplicates have been removed. If not, return a -1.
  if (newSize < minEigenVectorSize) {
    printf("    Chi-square test failed: eigenvalue array too small after duplicate removal. See the eigenvector output file.\n");
    return -1;
  }

  // We want to generate an average Chi-square value across various levels of
  // unfolding. Therefore, we iterate through the min and max unfolding pace
  // and then average the Chi-square values returned
  for (m = minUnfoldingPace; m < maxUnfoldingPace; m++) {

    // If the size / pace is fewer than 5 then skip this test
    if (newSize / m < 5) {
      continue;
    }

    // Get the Chi-Square value
    chiTest = getNNSDPaceChiSquare(newE, newSize, nnsdHistogramBin, m);
    avg_chiTest += chiTest;
    i++;

  }

  free(newE);

  // return the average Chi-square value
  return avg_chiTest / i;
}

/**
 * Performs a Chi-square test by comparing NNSD of the eigenvalues
 * to a Poisson distribution for a given pace.
 *
 * @param float* eigens
 *   The eigenvalue array
 * @param int size
 *   The size of the eigenvalue array
 * @param float bin
 *   The relative histogram bin size
 * @param int pace
 *   The unfolding pace
 * @param RMTParameters
 *
 *
 * @return float
 *   A Chi-square value, or -1 on failure
 */

float RMTThreshold::getNNSDPaceChiSquare(float* eigens,
    int size, float bin, int pace) {

  // The nearest neighbor spacing array.
  float * edif;
  float obj;
  float expect;
  float chi = 0;
  int i, j, count;


  // Unfolding will calculate the nearest neighbor spacing via estimation using
  // a spline curve. It returns a new array of length size - 1
  edif = unfolding(eigens, size, pace);
  size = size - 1;

  // Construct a histogram of (3.0/bin) + 1 bins.  If bin is 0.05 then the
  // number of bins is 60.  This corresponds to 60 degrees of freedom in the
  // Chi-square test. Therefore, if the desired p-value is 0.001 the
  // desired Chi-square value will be 99.607.  We don't actually save the
  // histogram, we just calculate the observed frequency and use that for
  // calculation of the Chi-square value.
  int n = (int) (3.0 / bin) + 1;
  for (i = 0; i < n; i++) {
    count = 0;

    // Create the histogram (or NNSD).  We are only interested in the
    // bins between 0 and 3.
    for (j = 0; j < size; j++) {
      if (edif[j] > i * bin && edif[j] < (i + 1) * bin) {
        count++;
      }
    }

    // https://books.google.com/books?id=Kp3Nx03_gMwC, pg 12.
    // The probability of s, p(s), from a poisson distribution is the integral
    // between a and b, where a and b are points between 0 and 3 (our
    // range of interest). The p(s) is calculated as:
    //
    //      p(s)ds = e(-s)ds
    //   or
    //      p(s) = e(-s) - e(-s+1)
    //
    //  The expected value is therefore the probability, p(s) times the
    //  number of eigenvalues.
    //
    expect = (exp(-1 * i * bin) - exp(-1 * (i + 1) *bin)) * size;

    // Perform the summation used for calculating the Chi-square value.
    // When the looping completes we will have the final Chi-square value.
    obj = (float) count;
    chi += (obj - expect) * (obj - expect) / expect;
  }
  free(edif);
  return chi;
}
