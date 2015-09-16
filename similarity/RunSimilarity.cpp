#include "RunSimilarity.h"

/**
 * Prints the command-line usage instructions for the similarity command
 */
void RunSimilarity::printUsage() {
  printf("\n");
  printf("Usage: ./kinc similarity [options]\n");
  printf("The list of required options:\n");
  printf("  --ematrix|-e      The tab delimited file name that contains the expression matrix.\n");
  printf("                    The rows must be genes or probe sets and columns are samples.\n");
  printf("                    If a header row is present it must only contain the list of\n");
  printf("                    genes (i.e. will be one column shorter than all other rows).\n");
  printf("  --rows|-r         The number of lines in the ematrix file including the header\n");
  printf("                    row if it exists\n");
  printf("  --cols|-c         The number of columns in the input file\n");
  printf("  --method|-m       The correlation method to use. Supported methods include\n");
  printf("                    Pearson's correlation ('pc'), Spearman's rank correlation ('sc')\n");
  printf("                    and Mutual Information ('mi'). Provide either 'pc', 'sc', or\n");
  printf("                    'mi' as values respectively.\n");
  printf("\n");
  printf("Optional Expression Matrix Arguments:\n");
  printf("  --omit_na         Provide this flag to ignore missing values.\n");
  printf("  --na_val|-n       A string representing the missing values in the input file\n");
  printf("                    (e.g. NA or 0.000)\n");
  printf("  --func|-f         A transformation function to apply to elements of the ematrix.\n");
  printf("                    Values include: log, log2 or log10. Default is to not perform\n");
  printf("                    any transformation.\n");
  printf("  --headers         Provide this flag if the first line of the matrix contains\n");
  printf("                    headers.\n");
  printf("\n");
  printf("Optional Similarity Arguments:\n");
  printf("  --min_obs|-o      The minimum number of observations (after missing values\n");
  printf("                    removed) that must be present to calculate a simililarity score.\n");
  printf("                    Default is 30.\n");
  printf("\n");
  printf("Optional Mutual Information Arguments:\n");
  printf("  --mi_bins|-b      Use only if the method is 'mi'. The number of bins for the\n");
  printf("                    B-spline estimator function for MI. Default is 10.\n");
  printf("  --mi_degree|-d  MT-ematrix-log-no.sc6.bin  Use only if the method is 'mi'. The degree of the\n");
  printf("                    B-spline estimator function for MI. Default is 3.\n");
  printf("\n");
  printf("Optional Clustering Arguments:\n");
  printf("  --clustering|-l   Perform pre-clustering of pair-wise samples prior to performing\n");
  printf("                    the similarity calculations. The similarity calculations will be\n");
  printf("                    performed on clusters greater than or equal to the min_obs\n");
  printf("                    argument. The value provided is the clustering method. The\n");
  printf("                    only method currently supported is mixture models \n");
  printf("                    (e.g. --clustering mixmod). Note: clustering significantly\n");
  printf("                    increases execution time and generates a different output format\n");
  printf("                    that can be hundres of gigabytes in size\n");
  printf("                    method to use. Currently only mixture models are supported.\n");
  printf("  --num_jobs|-j     Pair-wise sample clustering is highly parallel,\n");
  printf("                    thus the task can be split into multiple separate jobs. KINC \n");
  printf("                    can be run multiple times each working on a different portion \n");
  printf("                    of the total comparisions.  Use this argument to specify the \n");
  printf("                    total number of separate jobs that will be executed in parallel.\n");
  printf("                    KINC will not launch these. They must be launched manually, each\n");
  printf("                    with a different job_index. Defaults to 1\n");
  printf("  --job_index|-i    When the argument num_jobs is provided, this argument indicates which\n");
  printf("                    job is being run. For example, if num_jobs is set to 10, then this\n");
  printf("                    argument should be any number between 1 and 10. Defaults to 1.\n");
  printf("\n");
  printf("Optional Mixture Module Clustering Arguments:\n");
  printf("  --max_clusters|-a The maximum number of clusters that can be found for each\n");
  printf("                    pairwise comparision. Values between 2 to 10 are reasonable.\n");
  printf("  --criterion|-r    The Mixture module criterion to use. Valid values include:\n");
  printf("                    BIC, ICL, NEC, CV or DCV.  ICL may be more appropriate for\n");
  printf("                    smaller datasets (e.g. < 100 samples), and BIC for larger.\n");
  printf("\n");
  printf("For Help:\n");
  printf("  --help|-h       Print these usage instructions\n");
  printf("\n");
  printf("Note: similarity values are set to NaN if there weren't enough observations\n");
  printf("to perform the calculation.\n");
}
/**
 * The function to call when running the 'similarity' command.
 */
RunSimilarity::RunSimilarity(int argc, char *argv[]) {

  // Initialize some of the program parameters.
  clustering = NULL;
  min_obs = 30;

  // Initialize the clustering arguments.
  strcpy(criterion, "BIC");
  max_clusters = 5;
  num_jobs = 1;
  job_index = 0;

  // Defaults for mutual informatin B-spline estimate.
  mi_bins = 10;
  mi_degree = 3;

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
      {"method",       required_argument, 0,  'm' },
      {"min_obs",      required_argument, 0,  'o' },
      // Mutual information options.
      {"mi_bins",      required_argument, 0,  'b' },
      {"mi_degree",    required_argument, 0,  'd' },
      // Clustering options.
      {"clustering",   required_argument, 0,  'l' },
      {"num_jobs",     required_argument, 0,  'j' },
      {"job_index",    required_argument, 0,  'i' },
      // Mixture module options.
      {"criterion",    required_argument, 0,  't' },
      {"max_clusters", required_argument, 0,  'a' },
      // Expression matrix options.
      {"rows",         required_argument, 0,  'r' },
      {"cols",         required_argument, 0,  'c' },
      {"headers",      no_argument,       &headers,  1 },
      {"omit_na",      no_argument,       &omit_na,  1 },
      {"func",         required_argument, 0,  'f' },
      {"na_val",       required_argument, 0,  'n' },
      {"ematrix",      required_argument, 0,  'e' },
      // Last element required to be all zeros.
      {0, 0, 0,  0 }
    };

    // get the next option
    c = getopt_long(argc, argv, "m:o:b:d:j:i:t:a:l:r:c:f:n:e:h", long_options, &option_index);

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
      case 'o':
        min_obs = atoi(optarg);
        break;
      // Mutual information options.
      case 'b':
        mi_bins = atoi(optarg);
        break;
      case 'd':
        mi_degree = atoi(optarg);
        break;
      // Clustering options.
      case 'l':
        clustering = optarg;
        break;
      case 'j':
        num_jobs = atoi(optarg);
        break;
      case 'i':
        job_index = atoi(optarg);
        break;
      case 't':
        strcpy(criterion, optarg);
        break;
      case 'a':
        max_clusters = atoi(optarg);
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

  // Nake sure the similarity method is valid.
  if (!method) {
    fprintf(stderr,"Please provide the method (--method option).\n");
    exit(-1);
  }
  if (strcmp(method, "pc") != 0 &&
      strcmp(method, "mi") != 0 &&
      strcmp(method, "sc") != 0 ) {
    fprintf(stderr,"Error: The method (--method option) must either be 'pc', 'sc' or 'mi'.\n");
    exit(-1);
  }

  // make sure the required arguments are set and appropriate
  if (!infilename) {
    fprintf(stderr,"Please provide an expression matrix (--ematrix option).\n");
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
  // make sure the input file exists
  if (access(infilename, F_OK) == -1) {
    fprintf(stderr,"The input file does not exists or is not readable.\n");
    exit(-1);
  }

  // Validate the clustering options.
  if (clustering) {
    // Make sure the clustering method is valid.
    if (strcmp(clustering, "mixmod") != 0) {
      fprintf(stderr,"Error: The clustering type (--clustering option) must be 'mixmod'. It is currently the only clustering method supported.\n");
          exit(-1);
    }

    // If this is mixture model clustering the validate those options.
    if (strcmp(clustering, "mixmod") == 0) {
      // Make sure the MixMod criteria is correct.
      bool $mmc_is_good = false;
      if (strcmp(criterion, "BIC") == 0) {
        $mmc_is_good = true;
      }
      else if (strcmp(criterion, "ICL") == 0) {
        $mmc_is_good = true;
      }
      else if (strcmp(criterion, "NEC") == 0) {
        $mmc_is_good = true;
      }
      else if (strcmp(criterion, "CV") == 0) {
        $mmc_is_good = true;
      }
      else if (strcmp(criterion, "DCV") == 0) {
        $mmc_is_good = true;
      }
      if (!$mmc_is_good) {
        fprintf(stderr, "Error: The mixture model criterion must be one of: BIC, ICL, NEC, CV, DCV (--criterion option).\n");
        exit(-1);
      }
      // Make sure there is a valid max clusters value.
      if (max_clusters < 2 || max_clusters > 10 ){
        fprintf(stderr, "Error: Please select a maximum cluters between 2 and 10. (--max_clusters option).\n");
        exit(-1);
      }
      // Make sure the number of jobs is set
      if (num_jobs < 1) {
        fprintf(stderr, "Error: Please provide a positive integer for the number of jobs to run. (--num_jobs option).\n");
        exit(-1);
      }
      if (job_index < 1 || job_index > num_jobs) {
        fprintf(stderr, "Error: The job index must be between 1 and %d. (--job_index option).\n", num_jobs);
        exit(-1);
      }
    }
  }

  // Create and initialize the histogram for the distribution of coefficients.
  histogram = (int *) malloc(sizeof(int) * HIST_BINS + 1);
  for (int m = 0; m < HIST_BINS + 1; m++) {
    histogram[m] = 0;
  }

  // Retrieve the data from the EMatrix file.
  ematrix = new EMatrix(infilename, rows, cols, headers, omit_na, na_val, func);

  if (headers) {
    printf("  Skipping header lines\n");
  }
  printf("  Performing transformation: %s \n", func);
  if (omit_na) {
    printf("  Missing values are: '%s'\n", na_val);
  }
  printf("  Required observations: %d\n", min_obs);
  printf("  Using similarity method: '%s'\n", method);
  if (strcmp(method, "mi") ==0) {
    printf("  Bins for B-Spline estimate of MI: %d\n", mi_bins);
    printf("  Degree for B-Spline estimate of MI: %d\n", mi_degree);
  }
  if (clustering) {
    if (strcmp(clustering, "mixmod") == 0) {
      printf("\n  Mixture Model Settings:\n");
      printf("    Mixture model criterion: %s\n", criterion);
      printf("    Max clusters: %d\n", max_clusters);
      printf("    Number of jobs: %d\n", num_jobs);
      printf("    Job index: %d\n", job_index);
    }
  }
}

/**
 * Implements the destructor.
 */
RunSimilarity::~RunSimilarity() {
  delete ematrix;
  free(histogram);
}

/**
 *
 */
void RunSimilarity::execute() {
  if (clustering) {
    if (strcmp(clustering, "mixmod") == 0) {
      // TODO: the clustering and traditional looping should be more
      // integrated. It should be possible for MPI to support multiple
      // jobs for either method.
      MixtureModelClustering * mwc = new MixtureModelClustering(ematrix,
          min_obs, num_jobs, job_index, method, criterion, max_clusters);
      mwc->run();
    }
  }
  else {
    executeTraditional();
  }
}

/**
 *
 */
void RunSimilarity::executeTraditional() {
  // The output file name.
  char outfilename[1024];
  // Used for pairwise comparision of the same gene.
  float one = 1.0;
  // The number of binary files needed to store the matrix.
  int num_bins;
  // The current binary file number.
  int curr_bin;
  // Holds the number of rows in the file.
  int bin_rows;
  // The total number of pair-wise comparisions to be made.
  int total_comps;
  // The number of comparisions completed during looping.
  int n_comps;
  // The number of genes.
  int num_genes = ematrix->getNumGenes();
  // The binary output file prefix
  char * fileprefix = ematrix->getFilePrefix();
  char outdir[100];

  // calculate the number of binary files needed to store the similarity matrix
  num_bins = (num_genes - 1) / ROWS_PER_OUTPUT_FILE;

  // Make sure the output directory exists
  if (strcmp(method, "sc") == 0) {
    strcpy((char *)&outdir, "./Spearman");
  }
  if (strcmp(method, "pc") == 0) {
    strcpy((char *)&outdir, "./Pearson");
  }
  if (strcmp(method, "mi") == 0) {
    strcpy((char *)&outdir, "./MI");
  }
  struct stat st = {0};
  if (stat(outdir, &st) == -1) {
    mkdir(outdir, 0700);

  }

  total_comps = (num_genes * num_genes) / 2;
  n_comps = 0;

  // each iteration of m is a new output file
  printf("Calculating correlations...\n");
  for (curr_bin = 0; curr_bin <= num_bins; curr_bin++) {

    // calculate the limit on the rows to output based on where we are in the calculation
    if (curr_bin < num_bins) {
      bin_rows = (curr_bin + 1) * ROWS_PER_OUTPUT_FILE;
    }
    else {
      bin_rows = num_genes;
    }

    // the output file will be located in the Spearman directory and named based on the input file info
    sprintf(outfilename, "%s/%s.%s%d.bin", outdir,fileprefix, method, curr_bin);
    printf("Writing file %d of %d: %s... \n", curr_bin + 1, num_bins + 1, outfilename);
    FILE * outfile = fopen(outfilename, "wb");

    // write the size of the matrix
    fwrite(&num_genes, sizeof(num_genes), 1, outfile);
    // write the number of lines in this file
    int num_lines = bin_rows - (curr_bin * ROWS_PER_OUTPUT_FILE);
    fwrite(&num_lines, sizeof(num_lines), 1, outfile);

    // iterate through the genes that belong in this file
    for (int j = curr_bin * ROWS_PER_OUTPUT_FILE; j < bin_rows; j++) {
      // iterate through all the genes up to j (only need lower triangle)
      for (int k = 0; k <= j; k++) {
        n_comps++;
        if (n_comps % 1000 == 0) {
          statm_t * memory = memory_get_usage();
          printf("Percent complete: %.2f%%. Mem: %ldb. \r", (n_comps / (float) total_comps) * 100, memory->size);
          free(memory);
        }
        if (j == k) {
          // correlation of an element with itself is 1
          fwrite(&one, sizeof(one), 1, outfile);
          continue;
        }

        float score;
        PairWiseSet * pwset = new PairWiseSet(ematrix, j, k);

        // Perform the appropriate calculation based on the method
        if (strcmp(method, "pc") == 0) {
          PearsonSimilarity * pws = new PearsonSimilarity(pwset, min_obs);
          pws->run();
          score = (float) pws->getScore();
          delete pws;
        }
        else if(strcmp(method, "mi") == 0) {
          MISimilarity * pws = new MISimilarity(pwset, min_obs, mi_bins, mi_degree);
          pws->run();
          score = (float) pws->getScore();
          delete pws;
        }
        else if(strcmp(method, "sc") == 0) {
          SpearmanSimilarity * pws = new SpearmanSimilarity(pwset, min_obs);
          pws->run();
          score = (float) pws->getScore();
          delete pws;
        }

        delete pwset;
        fwrite(&score, sizeof(float), 1, outfile);

        // if the historgram is turned on then store the value in the correct bin
        if (score < 1 && score > -1) {
          if (score < 0) {
            score = -score;
          }
          histogram[(int)(score * HIST_BINS)]++;
        }
      }
    }
    fclose(outfile);
  }

  // Write the historgram
  writeHistogram();

  printf("\nDone.\n");
}

/**
 * Prints the histogram to a file.
 *
 * @param CCMParameters params
 *   An instance of the CCM parameters structsimilarity
 *
 * @param int * histogram
 *   An integer array created by the init_histogram function and
 *   populated by calculate_MI or calculate_person.
 *
 */
void RunSimilarity::writeHistogram() {
  char outfilename[1024];  // the output file name

  // output the correlation histogram
  sprintf(outfilename, "%s.%s.corrhist.txt", ematrix->getFilePrefix(), method);
  FILE * outfile = fopen(outfilename, "w");
  for (int m = 0; m < HIST_BINS; m++) {
    fprintf(outfile, "%lf\t%d\n", 1.0 * m / (HIST_BINS), histogram[m]);
  }
  fclose(outfile);
}

