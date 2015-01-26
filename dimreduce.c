#include "dimreduce.h"

int do_dimreduce(int argc, char *argv[], int mpi_id, int mpi_num_procs) {

  // variables used for timing of the software
  time_t start_time = time(0);
  time_t now;

  // the value returned by getopt_long
  int c;

  // The struct containing the input parameters
  static CCMParameters params;

  // initialize some of the program parameters
  params.perf = 1;
  params.omit_na = 0;
  params.headers = 0;
  params.rows = 0;
  params.cols = 0;
  params.do_log10 = 0;
  params.do_log2 = 0;
  params.do_log = 0;
  params.min_obs = 30;

  strcpy(params.func, "none");

  // loop through the incoming arguments until the
  // getopt_long function returns -1. Then we break out of the loop
  while(1) {
    int option_index = 0;

    // specify the long options. The values returned are specified to be the
    // short options which are then handled by the case statement below
    static struct option long_options[] = {
      {"omit_na",  no_argument,       &params.omit_na,  1 },
      {"perf",     no_argument,       &params.perf,     1 },
      {"headers",  no_argument,       &params.headers,  1 },
      {"help",     no_argument,       0,  'h' },
      {"ematrix",  required_argument, 0,  'e' },
      {"method",   required_argument, 0,  'm' },
      {"rows",     required_argument, 0,  'r' },
      {"cols",     required_argument, 0,  'c' },
      {"min_obs",  required_argument, 0,  'o' },
      {"func",     required_argument, 0,  'f' },
      {"na_val",   required_argument, 0,  'n' },
      {"mi_bins",  required_argument, 0,  'b' },
      {"mi_degree",required_argument, 0,  'd' },
      {0, 0, 0,  0 }  // last element required to be all zeros
    };
    // get the next option
    c = getopt_long(argc, argv, "e:r:c:m:o:n:f:h", long_options, &option_index);

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
      case 'r':
        params.rows = atoi(optarg);
        break;
      case 'c':
        params.cols = atoi(optarg);
        break;
      case 'o':
        params.min_obs = atoi(optarg);
        break;
      case 'n':
        params.na_val = optarg;
        break;
      case 'f':
        strcpy(params.func, optarg);
        break;
      case 'h':
        print_dimreduce_usage();
        exit(-1);
        break;
      case '?':
        exit(-1);
        break;
      case ':':
        print_dimreduce_usage();
        exit(-1);
        break;
      default:
        print_dimreduce_usage();
    }
  }

  // make sure the required arguments are set and appropriate
  if (!params.infilename) {
    fprintf(stderr, "Please provide an expression matrix (--ematrix option).\n");
    exit(-1);
  }
  // make sure we have a positive integer for the rows and columns of the matrix
  if (params.rows < 0 || params.rows == 0) {
    fprintf(stderr, "Please provide a positive integer value for the number of rows in the \n");
    fprintf(stderr, "expression matrix (--rows option).\n");
    exit(-1);
  }
  if (params.cols < 0 || params.cols == 0) {
    fprintf(stderr, "Please provide a positive integer value for the number of columns in\n");
    fprintf(stderr, "the expression matrix (--cols option).\n");
    exit(-1);
  }

  // make sure the input file exists
  if (access(params.infilename, F_OK) == -1) {
    fprintf(stderr, "Error: The input file does not exists or is not readable.\n");
    exit(-1);
  }

  if (params.omit_na && !params.na_val) {
    fprintf(stderr, "Error: The missing value string should be provided (--na_val option).\n");
    exit(-1);
  }

  if (params.headers) {
    printf("  Skipping header lines\n");
  }
  printf("  Performing transformation: %s \n", params.func);
  printf("  Required observations: %d\n", params.min_obs);
  if (params.omit_na) {
    printf("  Missing values are: '%s'\n", params.na_val);
  }

  // remove the path and extension from the filename
  char * filename = basename(params.infilename);
  strcpy(params.fileprefix, filename);
  char * p = rindex(params.fileprefix, '.');
  if (p) {
    p[0] = 0;
  }

  // if the function is log2 then set the do_log2 flag
  if (strcmp(params.func, "log10") == 0) {
    params.do_log10 = 1;
  }
  if (strcmp(params.func, "log2") == 0) {
    params.do_log2 = 1;
  }
  if (strcmp(params.func, "log") == 0) {
    params.do_log = 1;
  }

  // if we have a header line in the input file then
  // subtract one from the number of rows
  if (params.headers) {
    params.rows--;
  }

  EMatrix ematrix = load_ematrix(params);
  double ** data = ematrix.data;
  int i, j;

  // Open the clustering file for writing. Use the MPI ID in the filename.
  FILE ** fps = open_output_files(params, mpi_id);

  // Calculate the total number of comparisions and how many will be
  // performed by this process. We subtract 1 from the first params.rows
  // because we do not calculate the diagnoal.
  int num_rows = params.rows - 1;
  int total_comps  = num_rows * (num_rows + 1) / 2;
  int comps_per_process = total_comps / mpi_num_procs;
  int comp_start = mpi_id * comps_per_process;
  int comp_stop = mpi_id * comps_per_process + comps_per_process;

  // If this is the last process and there are some remainder comparisions
  // then we need to add them to the stop
  if (mpi_id + 1 == mpi_num_procs) {
    comp_stop = total_comps;
  }
  printf("%d. Performing %d comparisions\n", mpi_id + 1, comp_stop - comp_start);
  fflush(stdout);

  // Perform the pair-wise royston test and clustering
  int n_comps = 0;
  int my_comps = 0;
  for (i = 0; i < params.rows; i++) {
    for (j = 0; j < params.rows; j++) {

      // We only need to calculate royston in the lower triangle of the
      // full pair-wise matrix
      if (j >= i) {
        continue;
      }

      // If this computation is not meant for this process, then skip it.
      if (n_comps < comp_start || n_comps >= comp_stop) {
        n_comps++;
        continue;
      }

      if (n_comps % 100 == 0) {
        // Get the amount of memory used.
        statm_t * memory = memory_get_usage();

        // Calculate the number of days left.
        now = time(0);
        float seconds_passed = now - start_time;
        float minutes_passed = (seconds_passed) / 60.0;
        float percent_complete = (my_comps / (float) (comp_stop - comp_start)) * 100;
        float comps_per_minute = my_comps / minutes_passed;
        float total_time = (comp_stop - comp_start) / comps_per_minute;
        float minutes_left = total_time - minutes_passed;
        float hours_left = minutes_left / 60;
        float days_left = hours_left / 24;

        // Write progress report.
        printf("%d. Complete: %.4f%%. Mem: %ldb. Remaining: %.2fh; %.2fd. Coords: %d, %d.        \r",
          mpi_id + 1,
          percent_complete,
          memory->size,
          hours_left,
          days_left,
          i,
          j
        );
        free(memory);
      }
      n_comps++;
      my_comps++;

      double *a = data[i];
      double *b = data[j];

      // Initialize the arrays that will be used for containing a and b
      // but with missing values removed.
      double * a2 = (double *) malloc(sizeof(double) * params.cols);
      double * b2 = (double *) malloc(sizeof(double) * params.cols);

      // The kept variable is an array of zeros and ones and is the size of the
      // number of samples.  A zero indicates the sample is not included in
      // the pair-wise comparision and a one indicates it is.
      int * kept = (int *) malloc(sizeof(int) * params.cols);
      int n2;

      // Remove any missing values before calculating Royston's H test.
      // This call will return a2 and b2 which are the same as a and b but
      // with missing values removed. n2 gets set to the size of a2 and b2.
      remove_missing_paired(a, b, params.cols, a2, b2, &n2, kept);

      // Perform the clustering if we have enough samples.
      if (n2 > 0) {

        // Perform the clustering.
        PairWiseClusters * clusters = clustering(a2, i, b2, j, n2, ematrix, params, 0.075, 0);

        // update the clusters to include zeros for any samples with a missing value.
        update_pairwise_cluster_samples(kept, params.cols, clusters);

        // Write the clusters to the file.
        write_pairwise_cluster_samples(clusters, fps);

        // clean up the memory.
        free_pairwise_cluster_list(clusters);
      }
      // One of these genes has no observations so, just write out an
      // empty cluster.
      else {
        // Create an empty clusters set
        PairWiseClusters * new = new_pairwise_cluster_list();
        new->gene1 = i;
        new->gene2 = j;
        new->num_samples = params.cols;
        new->samples = kept;
        new->next = NULL;
        new->cluster_size = 0;
        new->pcc = NAN;
        write_pairwise_cluster_samples(new, fps);
        free(new);
      }

      // Release the memory for a2 and b2.
      free(a2);
      free(b2);
      free(kept);
    }
  }

  close_output_files(fps);

  return 1;
}

/**
 * @param double * a2
 *   A row from the ematrix with missing values removed
 * @param int x
 *   The index of a2 in the ematrix
 * @param double * b2
 *   A row from the ematrix with missing values removed.
 * @param in y
 *   The index of b2 in the ematrix
 * @param int n2
 *   The size of both a2 and b2.
 * @param Ematrix ematrix
 *   The two-dimensional ematrix where rows are genes and columns are samples.
 * @param CCMParameters params
 *   The parameters provided to the program
 * @param float bw
 *   The bandwith argument for mean shift clustering.  Default should be 0.075.
 * @param int level
 *   An integer indicating the recursion level. It should always be set to
 *   zero by the caller of this function.
 */
PairWiseClusters * clustering(double *a2, int x, double *b2, int y, int n2,
    EMatrix ematrix, CCMParameters params, float bw, int level) {

  // Variables used for looping.
  int k, l, j;
  int nkept;
  int * ckept;

  PairWiseClusters * result = new_pairwise_cluster_list();

  // Perform mean shift clustering (MSC)
  MeanShiftClusters clusters;
  clusters = meanshift2D(a2, b2, n2, bw);

  // Count the number of clusters that are larger than min_obs
  int num_large = 0;
  for(k = 0; k < clusters.num_clusters; k++) {
    if (clusters.sizes[k] >= params.min_obs) {
      num_large++;
    }
  }

  // Iterate through all of the clusters.
  for(k = 0; k < clusters.num_clusters; k++) {

    // For easier access create variables for the current cluster size.
    int size = clusters.sizes[k];

    // Create separate vectors for the x and y coordinates.
    double *cx = (double *) malloc(sizeof(double) * size);
    double *cy = (double *) malloc(sizeof(double) * size);
    for (j = 0; j < size; j++) {
      cx[j] = 0;
      cy[j] = 0;
    }
    // Add the values to the cx & cy arrays
    int l = 0;
    for (j = 0; j < n2; j++) {
      if (clusters.cluster_label[j] == k + 1) {
        cx[l] = a2[j];
        cy[l] = b2[j];
        l++;
      }
    }

    // Discover any outliers for clusters with size >= min_obs
    Outliers outliersCx;
    Outliers outliersCy;
    if (clusters.sizes[k] >= params.min_obs) {
      outliersCx = outliers_iqr(cx, size, 1.5);
      outliersCy = outliers_iqr(cy, size, 1.5);
    }
    else {
      outliersCx.n = 0;
      outliersCy.n = 0;
    }

    // Create an array of the kept samples for this cluster. Don't include
    // any samples whose coordinates are considered outliers or who are not
    // in the cluster.
    ckept = (int *) malloc(sizeof(int) * n2);
    for (l = 0; l < n2; l++) {
      ckept[l] = 0;
    }
    nkept = 0;
    for (l = 0; l < n2; l++) {
      // Is this sample is in the cluster? if so then also make sure it's
      // not an outlier.
      if (clusters.cluster_label[l] == k + 1) {
        // Iterate through the outlier points and compare to this
        // sapmle's points.  If there is a match then mark as an outlier
        // and exclude it from the samples that are kept.  First check the
        // x coordinate
        int is_outlier = 0;
        for (j = 0; j < outliersCx.n; j++) {
          if (a2[l] == outliersCx.outliers[j]) {
            is_outlier = 1;
            break;
          }
        }
        // Second check the y coordinate.
        for (j = 0; j < outliersCy.n; j++) {
          if (b2[l] == outliersCy.outliers[j]) {
            is_outlier = 1;
            break;
          }
        }
        // if it's not an outlier then keep it.
        if (!is_outlier) {
          ckept[l] = 1;
          cx[nkept] = a2[l];
          cy[nkept] = b2[l];
          nkept++;
        }
      }
    }
    if (outliersCx.n > 0) {
      free(outliersCx.outliers);
    }
    if (outliersCy.n > 0) {
      free(outliersCy.outliers);
    }

    // Now after we have removed outliers and non-cluster samples,
    // makes sure we still have the minimum observations.
    if (nkept >= params.min_obs) {

      // If the recursion level is greater than zero the we've already clustered
      // at least once.  When we reach this point we've clustered again. If
      // we only have a single cluster at this point then we can't cluster
      // anymore and we can perform Spearman.
      if (level > 0 && num_large == 1) {
        // Initialize the workspace needed for Spearman's calculation.
        double workspace[2 * params.rows];
        float rho = NAN;
        rho = gsl_stats_spearman(cx, 1, cy, 1, nkept, workspace);
        PairWiseClusters * new = new_pairwise_cluster_list();
        new->gene1 = x;
        new->gene2 = y;
        new->num_samples = n2;
        new->samples = ckept;
        new->next = NULL;
        new->cluster_size = nkept;
        new->pcc = rho;
        add_pairwise_cluster_list(result, new);
      }

      // If this is the first level of clsutering then we want to cluster
      // again, but only if there is one or more large cluster.
      // We do this because we haven't yet settled on a distinct
      // "expression mode".
      else {
        PairWiseClusters * children = clustering(cx, x, cy, y, nkept, ematrix, params, 0.09, level + 1);
        update_pairwise_cluster_samples(ckept, n2, children);
        add_pairwise_cluster_list(result, children);
      }
    }
    // If the cluster is too small then just add it without
    // correlation analysis or further sub clustering.
    else {
      PairWiseClusters * new = new_pairwise_cluster_list();
      new->gene1 = x;
      new->gene2 = y;
      new->num_samples = n2;
      new->samples = ckept;
      new->next = NULL;
      new->cluster_size = nkept;
      new->pcc = NAN;
      add_pairwise_cluster_list(result, new);
    }
    free(cx);
    free(cy);
  }


  // free the memory
  free(clusters.cluster_label);
  for(k = 0; k < clusters.num_clusters; k++) {
    for (l = 0; l < clusters.sizes[k]; l++) {
      free(clusters.clusters[k][l]);
    }
    free(clusters.clusters[k]);
  }
  free(clusters.sizes);
  free(clusters.clusters);

  return result;
}


/**
 * Intializes the head PairWiseCluster object.
 */
PairWiseClusters * new_pairwise_cluster_list() {
  PairWiseClusters * pws = (PairWiseClusters *) malloc(sizeof(PairWiseClusters));
  pws->gene1 = -1;
  pws->gene2 = -1;
  pws->next = NULL;
  pws->num_samples = 0;
  pws->cluster_size = 0;
  pws->pcc = 0;

  return pws;
}

/**
 * Frees up all the memory in a PairWiseClusters object list.
 */
void free_pairwise_cluster_list(PairWiseClusters * head) {

  PairWiseClusters * curr = (PairWiseClusters *) head;
  PairWiseClusters * next = (PairWiseClusters *) head->next;
  if (curr->num_samples > 0) {
    free(curr->samples);
  }
  free(curr);
  while (next != NULL) {
    curr = next;
    next = (PairWiseClusters *) next->next;
    if (curr->num_samples > 0) {
      free(curr->samples);
    }
    free(curr);
  }
}
/**
 * Adds a new PairWiseClusters object to the list
 *
 * @param PairWiseClusters head
 *   The head object in the list
 * @param PairWiseClusters new
 *   The new object to add to the end of the list
 */
void add_pairwise_cluster_list(PairWiseClusters *head, PairWiseClusters *new) {

  // Check the list to see if it is empty. If so, then make this item the
  // new head.
  if (head->gene1 == -1) {
    head->gene1 = new->gene1;
    head->gene2 = new->gene2;
    head->next = new->next;
    head->samples = new->samples;
    head->num_samples = new->num_samples;
    head->cluster_size = new->cluster_size;
    head->pcc = new->pcc;
    return;
  }

  // Traverse the list to the end and then add the new item
  PairWiseClusters * curr;
  curr = head;
  while (curr->next != NULL) {
    curr = (PairWiseClusters * ) curr->next;
  }
  curr->next = (struct PairWiseClusters *) new;
}
/**
 * Updates the samples vector of a PairWiseClusters object.
 *
 * Because the clustering() function is recursive it is successivly called
 * with smaller and smaller sample sets.  When it returns it provdes  PWC
 * object with a samples array.  In the samples array, samples that are present
 * in the cluster are marked with a 1 and those not in the cluster are
 * set to 0.  The order of values correspondes to the sample order. The samples
 * array, however, only contains 1's and 0's for the samples provided to the
 * clustering() function. Therefore, the results need to be merged back into
 * the larger samples array.  This function does that.
 *
 * @param int *parent_samples
 *   The list of parent samples. It contains a list of 0's and 1's indicating
 *   which samples are to be kept.  Any samples not also found in the
 *   'new.samples' argument are set to 0.
 * @param int n
 *   The size of the parent_samples array.
 * @param new
 *   The new PWC object.
 */
void update_pairwise_cluster_samples(int * parent_samples, int n, PairWiseClusters * head) {

  PairWiseClusters * curr = (PairWiseClusters *) head;
  PairWiseClusters * next = (PairWiseClusters *) head->next;

  while (curr != NULL) {
    int z;
    int w = 0;

    // Create a new kept array that will update the current samples array.
    int * new_samples = (int *) malloc(sizeof(int) * n);

    // Iterate through all of the elements of the parent samples list.
    for (z = 0; z < n; z++) {
      // If the element is 1, meaning the sample is present int he parent
      // cluster, then check the current sample to see if was preserved during
      // sub clustering.
      if (parent_samples[z] == 1) {
        if (curr->samples[w] == 0) {
          new_samples[z] = 0;
        }
        // The element is kept in the pwc so preserve the 1 it in the new array.
        else {
          new_samples[z] = 1;
        }
        w++;
      }
      // The element is not kept originally so preserve the 0 in the new array.
      else {
        new_samples[z] = 0;
      }
    }
    // Free up the old samples array and replace it with a new one.
    free(curr->samples);
    curr->samples = new_samples;
    curr->num_samples = n;

    curr = next;
    if (next != NULL) {
      next = (PairWiseClusters *) next->next;
    }
  }
}

/**
 * Adds a line to the clustering file.
 *
 * The clustering file is used by KINC during pair-wise correlation analysis
 * to restrict which samples are used.  The file is tab delimited.
 * The format of the file is tab delimited with the following columns:
 *
 *   1)  gene 1 name
 *   2)  gene 2 name
 *   3)  cluster name.  A 0 indicates no clustering was performed.
 *   4)  a string of 0 and 1s indicating which samples to include when
 *       performing pair-wise comparisons.
 *
 * @param PairWiseSets pws
 *   The details for this line
 * @param FILE *cf
 *   The file pointer of the clustering file.
 */
void write_pairwise_cluster_samples(PairWiseClusters * pwc, FILE ** fps) {

  // The file pointer of the file to write to.
  FILE *fp;

  // Do nothing if the object is empty.
  if (pwc->gene1 == -1) {
    return;
  }

  // Iterate through the list of clusters and print each one.
  PairWiseClusters * curr = pwc;
  int cluster_id = 0;
  while (curr != NULL) {
    // Determine which file to write the output into
    if (isnan(curr->pcc)) {
      fp = fps[101];
    }
    else {
      float i1 = curr->pcc * 100.0;
      float i2 = fabs(i1);
      int i3 = (int) i2;
      fp = fps[i3];
    }
    if (curr->gene1 != -1) {
      fprintf(fp, "%i\t%i\t%i\t%i\t", curr->gene1, curr->gene2, curr->cluster_size, cluster_id + 1);
      int i;
      for (i = 0; i < curr->num_samples; i++) {
        fprintf(fp, "%i", curr->samples[i]);
      }
      fprintf(fp, "\t%f", curr->pcc);
      fprintf(fp, "\n");
      cluster_id++;
    }
    curr = (PairWiseClusters *) curr->next;
    fflush(fp);
  }
}

/**
 * Opens and creates 102 files for storing the clusters with correlation values.
 * Each file stores a range of 1/100 Spearman correlation values.
 */
FILE ** open_output_files(CCMParameters params, int mpi_id) {

  // make sure the output directory exists
  struct stat st = {0};
  if (stat("./clusters", &st) == -1) {
      mkdir("./clusters", 0700);
  }

  // Open up 102 files, one for a range of spearman correlation values
  FILE ** fps = malloc(sizeof(FILE *) * 102);

  int i =  0;
  char filename[1025];
  for (i = 0; i <= 100; i++) {
    sprintf(filename, "./clusters/%s.clusters.%03d.%03d.txt", params.fileprefix, i, mpi_id + 1);
    fps[i] = fopen(filename, "w");
  }

  sprintf(filename, "./clusters/%s.clusters.nan.%03d.txt", params.fileprefix, mpi_id + 1);
  fps[i] = fopen(filename, "w");

  return fps;
}

/**
 * Closes the 102 files that were opened.
 */
void close_output_files(FILE** fps) {
  int i =  0;
  for (i = 0; i <= 101; i++) {
    FILE * fp = fps[i];
    fprintf(fp, "#Done\n");
    fclose(fp);
  }
}

/**
 * Prints the command-line usage instructions for the similarity command
 */
void print_dimreduce_usage() {
  printf("\n");
  printf("Usage: ./kinc dimreduce [options]\n");
  printf("The list of required options:\n");
  printf("  --ematrix|-e The file name that contains the expression matrix.\n");
  printf("                 The rows must be genes or probe sets and columns are samples\n");
  printf("  --rows|-r    The number of lines in the input file including the header\n");
  printf("                 column if it exists\n");
  printf("  --cols|-c    The number of columns in the input file minus the first\n");
  printf("                 column that contains gene names\n");
  printf("\n");
  printf("Optional:\n");
  printf("  --omit_na     Provide this flag to ignore missing values. Defaults to 0.\n");
  printf("  --na_val|-n   A string representing the missing values in the input file\n");
  printf("                  (e.g. NA or 0.000)\n");
  printf("  --min_obs|-o  The minimum number of observations (after missing values\n");
  printf("                  removed) that must be present to perform correlation.\n");
  printf("                  Default is 30.\n");
  printf("  --func|-f     A transformation function to apply to elements of the ematrix.\n");
  printf("                  Values include: log, log2 or log10. Default is to not perform\n");
  printf("                  any transformation.\n");
  printf("  --perf        Provide this flag to enable performance monitoring.\n");
  printf("  --headers     Provide this flag if the first line of the matrix contains\n");
  printf("                  headers.\n");
  printf("  --help|-h     Print these usage instructions\n");
  printf("\n");
}
