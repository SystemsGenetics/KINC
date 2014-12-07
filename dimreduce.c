#include "dimreduce.h"

int do_dimreduce(int argc, char *argv[]) {

  //time_t start_time, end_time;  // variables used for timing of the software
  int c;                        // the value returned by getopt_long

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
    fprintf(stderr,"Please provide an expression matrix (--ematrix option).\n");
    exit(-1);
  }
  // make sure we have a positive integer for the rows and columns of the matrix
  if (params.rows < 0 || params.rows == 0) {
    fprintf(stderr,"Please provide a positive integer value for the number of rows in the \n");
    fprintf(stderr,"expression matrix (--rows option).\n");
    exit(-1);
  }
  if (params.cols < 0 || params.cols == 0) {
    fprintf(stderr,"Please provide a positive integer value for the number of columns in\n");
    fprintf(stderr,"the expression matrix (--cols option).\n");
    exit(-1);
  }

  // make sure the input file exists
  if (access(params.infilename, F_OK) == -1) {
    fprintf(stderr,"Error: The input file does not exists or is not readable.\n");
    exit(-1);
  }

  if (params.omit_na && !params.na_val) {
    fprintf(stderr,"Error: The missing value string should be provided (--na_val option).\n");
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
  int total_comps;  // the total number of pair-wise comparisons to be made
  int n_comps;      // the number of comparisons completed during looping
  int i, j;

  // Open the clustering file for writing.
  sprintf(filename, "%s.clusters.txt", params.fileprefix);
  FILE * cf = fopen(filename, "w");

  // Perform the pair-wise royston test and clustering
  total_comps = (params.rows * params.rows) / 2;
  for (i = 0; i < params.rows; i++) {
    for (j = 0; j < params.rows; j++) {

      n_comps++;
      if (n_comps % 1000 == 0) {
        printf("Percent complete: %.2f%%\r", (n_comps / (float) total_comps) * 100);
      }

      // We only need to calculate royston in the lower triangle of the
      // full pair-wise matrix
      if (j >= i) {
        continue;
      }
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
      if (n2 > params.min_obs) {
        // Perform the clustering
        PairWiseClusters * pws = clustering(a2, i, b2, j, n2, ematrix, params);

        // So, we have a good cluster. But, the samples in the pws object
        // only represent data from the subset of samples that formed the
        // cx and cy vectors. We need to adjust the samples array accordingly
        // so that it excludes anthing but the samples in the returned cluster.
        update_pairwise_cluster_samples(kept, params.cols, pws);

        // Now write out any clusters
        write_pairwise_cluster_samples(pws, cf);
      }

      // Release the memory for a2 and b2.
      free(a2);
      free(b2);
      free(kept);
    }
  }

  // Close the clustering file.
  fclose(cf);
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
 */
PairWiseClusters * clustering(double *a2, int x, double *b2, int y, int n2,
    EMatrix ematrix, CCMParameters params) {

  // Initialize the PWC struct to be empty (has -1 for genes and NULL for next)
  PairWiseClusters * result = new_pairiwse_cluster_list();

  // Variables used for looping.
  int k, l, j;
  int nkept;
  int * ckept = (int *) malloc(sizeof(int) * ematrix.num_samples);

  // Calculate Royston's H test for multivariate normality.
  double pv = royston2D(a2, b2, n2);
  //printf("(%d, %d), pv: %e\n", x + 1, y + 1, pv);

  // If the Royston's H test has p-value <= 0.05 it means
  // it does not appear to be multivariate normal, then do mean shift
  // clustering to cluster the measured points.
  if (pv > 0.05) {
    // Since this dataset is multivariate normal, return all samples
    // as being kept (set all elements of ckept to 1).
    for (l = 0; l < ematrix.num_samples; l++) {
      ckept[l] = 1;
    }

    // Return the entire sample set as a good cluster.
    result.gene1 = x;
    result.gene2 = y;
    result.num_samples = ematrix.num_samples;
    result.samples = ckept;
    result.cluster_label = 0;  // 0 indicates no clustering was performed.
    result.next = NULL;

    return result;
    //write_reduced_ematrix_line(pws, cf);
  }

  // Calculate the clusters.
  MeanShiftClusters clusters;
  clusters = meanshift2D(a2, b2, n2, 0.075);

  // Is there a cluster with the minimum observations?
  for(k = 0; k < clusters.num_clusters; k++) {

    // Skip clusters that are too small.
    if (clusters.sizes[k] < params.min_obs) {
      continue;
    }

    // For easier access create variables for the current cluster.
    double ** cluster = clusters.clusters[k];
    int size = clusters.sizes[k];

    // Create seperate vectors for the x and y coordinates. We
    // will look for outliers in both vectors
    double *cx = (double *) malloc(sizeof(double) * size);
    double *cy = (double *) malloc(sizeof(double) * size);
    for(j = 0; j < size; j++) {
      cx[j] = cluster[j][0];
      cy[j] = cluster[j][1];
    }

    // Discover any outliers in this cluster.
    Outliers outliersCx = outliers_iqr(cx, size, 1.5);
    Outliers outliersCy = outliers_iqr(cy, size, 1.5);

    // Re-initailize the cx and cy vectors to hold the new cluster (with
    // outliers removed)
    for(j = 0; j < size; j++) {
      cx[j] = 0;
      cy[j] = 0;
    }

    // Initialize the ckept array.
    for (l = 0; l < ematrix.num_samples; l++) {
      ckept[l] = 0;
    }

    // Create an array of the kept samples for this cluster. Don't include
    // any samples whose coordinates are considered outliers or who are not
    // in the cluster.
    nkept = 0;
    for (l = 0; l < ematrix.num_samples; l++) {
      // Is this sample is in the cluster, if so then also make sure it's
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

    free(outliersCx.outliers);
    free(outliersCy.outliers);

    // Now after we have removed outliers and non-cluster samples,
    // makes sure add_pairwise_cluster_to_listwe still have the minimum observations
    if (nkept >= params.min_obs) {
      PairWiseClusters pws;

      // Before we write out this cluster to the cluster file, check to make
      // sure that this cluster is multivariate normal. To do this, recursively
      // call this function.  This cluster will be cotinually clustered itself
      // until it is either multivariate normal or too small to cluster further.
      // if no valid cluster is found then the genes are set at -1.
      pws = clustering(cx, x, cy, y, nkept, ematrix, params);
      if (pws.gene1 > -1) {
        // So, we have a good cluster. But, the samples in the pws object
        // only represent data from the subset of samples that formed the
        // cx and cy vectors. We need to adjust the samples array accordingly
        // so that it excludes anthing but the samples in the returned cluster.
        update_pairwise_cluster_samples(ckept, nkept, pws);

        // Add the resulting cluster to the list
        add_pairiwse_cluster_list(&result, &pws);
      }
    }
    free(ckept);
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

  // Return the final list
  return result;
}

/**
 * Intializes the head PairWiseCluster object.
 */
PairWiseClusters * new_pairiwse_cluster_list() {
  PairWiseClusters * pws = (PairWiseClusters *) malloc(sizeof(PairWiseClusters));
  pws->gene1 = -1;
  pws->gene2 = -1;
  pws->next = NULL;

  return pws;
}
/**
 * Adds a new PairWiseClusters object to the list
 *
 * @param PairWiseClusters head
 *   The head object in the list
 * @param PairWiseClusters new
 *   The new object to add to the end of the list
 */
void add_pairiwse_cluster_list(PairWiseClusters *head, PairWiseClusters *new) {
  // If this is the first time we've added to the list then set this new
  // object to be the head.
  if (head->gene1 == -1) {
    head = new;
    return;
  }

  // Traverse the list to the end and then add the new item
  PairWiseClusters * temp;
  temp = head;
  while (temp->next != NULL) {
    temp = temp->next;
  }
  temp->next = new;
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
 * @param int *ckept
 *   The larger samples array. It contains a list of 0's and 1's indicating
 *   which samples are to be kept.  Any samples not also found in the
 *   'new.samples' argument are set to 0.
 * @param int nkept
 *   The size of the ckept array.
 * @param new
 *   The new PWC object.
 */
void update_pairwise_cluster_samples(int * ckept, int nkept, PairWiseClusters * new) {
  int z;
  int w = 0;
  for (z = 0; z < nkept; z++) {
    if (ckept[z] == 1) {
      if (new->samples[w] == 0) {
        ckept[z] = 0;
      }
      w++;
    }
  }
  free(new->samples);
  new->samples = ckept;
  new->num_samples = nkept;
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
void write_pairwise_cluster_samples(PairWiseClusters * pwc, FILE * cf) {

  // Do nothing if the object is empty.
  if (pwc.gene1 == -1) {
    return;
  }

  // Iterate through the list of clusters and print each one.
  PairWiseClusters * temp = pwc;
  while (temp->next != NULL) {
    fprintf(cf, "%i\t%i\t%i\t", pwc->gene1, pwc->gene2, pwc->cluster_label);
    int i;
    for (i = 0; i < pwc.num_samples; i++) {
      fprintf(cf, "%i", pwc.samples[i]);
    }
    fprintf(cf, "\n");
    temp = temp->next;
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
