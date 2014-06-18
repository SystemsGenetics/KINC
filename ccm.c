/** 
 * Description:
 * ------------
 * Calculates a pearson correlation coefficient matrix from an n x m 
 * expression matrix where the columns are the samples and rows are
 * the genes or probesets. Cells represent expression levels. The first
 * column of the expression matrix should be the name of the gene or
 * probeset.
 *
 *
 * Change Log:
 * ----------
 * 10/2014
 * 1) Input files no longer require an implied .txt extesion but the full filename
 *    should be specified
 * 2) Removed the spaces in the output file name and renamed file.
 * 3) Incorporated GSL Pearson calculation function because the pre-calculating
 *    of sum of row and sum of squared of the row couldn't account for missing values
 *
 * History:
 * --------
 * 04/2014
 * Stephen Ficklin added support for missing values and added comments to the code
 * and changed it so the '.txt' extension is no longer implied for incoming
 * file names
 *
 * 10/2011 
 * Created by Scott Gibson at Clemson University under direction of 
 * Dr. Melissa Smith in Collaboration with Alex Feltus, Feng Luo and Stephen
 * Ficklin
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <libgen.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <getopt.h>
#include <gsl/gsl_statistics.h>
#include "MutualInformation.h"

/**
 * Function Prototypes
 */
double ** load_ematrix();
int * init_histogram();
void print_histogram(int * histogram);
int calculate_MI(double ** data, int * histogram);
int calculate_pearson(double ** data, int * histogram);
void print_usage();

/**
 * Constants
 */

// a global variable for the number of rows in each output file
static int ROWS_PER_OUTPUT_FILE = 10000;

// the number of bins in the correlation value historgram
static int HIST_BINS = 100;


/**
 * Globals
 */

// Command-line option flags
static int perf;          // indicates if performance monitoring should be enabled
static int omit_na;       // indicates if missing values should be ignored
static int headers;

// Command-line option values
static int rows = 0;           // the number of rows in the expression matrix
static int cols = 0;           // the number of columns in the expression matrix
static char *infilename = 0;   // the input file name
static char *na_val = 0;       // specifies the value that represents a missing value
static char func[10] = "none"; // specifies the transformation function: log2, none
static char method[10] = "cor";// specifies the method: cor, mi
static int min_obs = 30;       // the minimum number of observations to calculate correlation


// other global variables
static int do_log10 = 0;      // set to 1 to perform log10 transformation
static int do_log2 = 0;       // set to 1 to perform log2 transformation
static int do_log  = 0;       // set to 1 to perform log transformation
static char fileprefix[1024]; // the input filename without the prefix

/**
 * The main subroutine.  Parses the input parameters and executes the program
 * accordingly.
 */
int main(int argc, char *argv[]) {

  time_t start_time, end_time;  // variables used for timing of the software
  int c;                        // the value returned by getopt_long
  int option_index = 0;         // the option index set by getopt_long

  // loop through the incoming arguments until the
  // getopt_long function returns -1. Then we break out of the loop
  while(1) {
    int option_index = 0;

    // specify the long options. The values returned are specified to be the
    // short options which are then handled by the case statement below
    static struct option long_options[] = {
      {"help",    no_argument,       0,  'h' },
      {"omit_na", no_argument,       &omit_na,  1 },
      {"perf",    no_argument,       &perf,     1 },
      {"headers", no_argument,       &headers,  1 },
      {"ematrix", required_argument, 0,  'e' },
      {"method",  required_argument, 0,  'm' },
      {"rows",    required_argument, 0,  'r' },
      {"cols",    required_argument, 0,  'c' },
      {"min_obs", required_argument, 0,  'o' },
      {"func",    required_argument, 0,  'f' },
      {"na_val",  required_argument, 0,  'n' },
      {0, 0, 0,  0 }  // last element required to be all zeros
    };

    // get the next option
    c = getopt_long(argc, argv, "he:r:c:m:o:n:f:", long_options, &option_index);

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
        infilename = optarg;
        break;
      case 'r':
        rows = atoi(optarg);
        break;
      case 'c':
        cols = atoi(optarg);
        break;
      case 'm':
        strcpy(method, optarg);
        printf("  Using method: '%s'\n", method);
        break;
      case 'o':
        min_obs = atoi(optarg);
        printf("  Requiring %d observations to perform pairwise calculations\n", min_obs);
        break;
      case 'n':
        na_val = optarg;
        printf("  Excluding missing values marked as '%s'\n", na_val);
        break;
      case 'f':
        strcpy(func, optarg);
        printf("  Performing %s transformation on expression matrix\n", func);
        break;
      case 'h':
        print_usage();
        exit(0);
      case '?':
        exit(-1);
        break;
      case ':':
        print_usage();
        exit(-1);
        break;
      default:
        print_usage();
    }
  }

  if (omit_na) {
    printf("  Excluding missing values\n");
  }
  if (headers) {
    printf("  Skipping header\n");
  }

  // make sure the required arguments are set and appropriate
  if (!infilename) {
    printf("Please provide an expression matrix (--ematrix option). Use the -h option for help.\n");
    exit(-1);
  }

  if (!method) {
    printf("Please provide the method (--method option). Use the -h option for help.\n");
    exit(-1);
  }

  // make sure we have a positive integer for the rows and columns of the matrix
  if (rows < 0 || rows == 0) {
    printf("Please provide a positive integer value for the number of rows in the \n");
    printf("expression matrix (--rows option). Use the -h option for help.\n");
    exit(-1);
  }
  if (cols < 0 || cols == 0) {
    printf("Please provide a positive integer value for the number of columns in\n");
    printf("the expression matrix (--cols option). Use the -h option for help.\n");
    exit(-1);
  }

  // make sure the input file exists
  if (access(infilename, F_OK) == -1) {
    printf("Error: The input file does not exists or is not readable.\n");
    exit(-1);
  }

  if (strcmp(method, "pc") != 0 && strcmp(method, "mi") != 0) {
    printf("Error: The method (--method option) must either be 'pc' or 'mi'. Use the -h option for help.\n");
    exit(-1);
  }

  // remove the path and extension from the filename
  char * filename = basename(infilename);
  strcpy(fileprefix, filename);
  char * p = rindex(fileprefix, '.');
  if (p) {
    p[0] = 0;
  }

  // if the function is log2 then set the do_log2 flag
  if (strcmp(func, "log10") == 0) {
    do_log10 = 1;
  }
  if (strcmp(func, "log2") == 0) {
    do_log2 = 1;
  }
  if (strcmp(func, "log") == 0) {
    do_log = 1;
  }

  // if we have a header line in the input file then subtract one from the number of rows
  if (headers) {
    rows--;
  }

  // if output of the histogram is enabled then initialize it
  int * histogram;
  histogram = init_histogram();

  // if performance monitoring is enabled the set the start time
  if (perf) {
    time(&start_time);
  }

  double ** data = load_ematrix();

  if (strcmp(method, "pc") == 0) {
    calculate_pearson(data, histogram);
  }
  else if(strcmp(method, "mi") == 0) {
    calculate_MI(data, histogram);
  }

  // if performance monitoring is enabled then write the timing data
  if(perf) {
    time(&end_time);
    FILE * timingfile = fopen("timingdata.txt", "a");
    fprintf(timingfile, "CCM Runtime with %d x %d %s input dataset: %.2lf min\n", rows, cols, infilename, difftime(end_time, start_time)/60.0);
  }

  // if output of the histogram is enabled then print it
  print_histogram(histogram);

  printf("Done.\n");
}

/**
 *
 */
double ** load_ematrix() {

  FILE *infile; // pointers to the input
  char gene_name[255];   // holds the gene or probe set name from the first column
  double ** data;
  int i, j;  // integers for looping

  // allocate the data array for storing the input expression matrix
  data = (double**) malloc(sizeof(double *) * rows);
  for (i = 0; i < rows; i++) {
    data[i] = malloc(sizeof(double) * cols);
  }

  // iterate through the lines of the expression matrix
  printf("Reading input file '%s'...\n", infilename);
  infile = fopen(infilename, "r");
  for (i = 0; i < rows; i++) {

    // skip over the header if one is provided
    if (i == 0 && headers) {
      printf("Skipping headers...\n");
      char element[1024];
      for (j = 0; j < cols; j++) {
        fscanf(infile, "%s", element);
      }
    }

    // the first entry on every line is a label string - read that in before the numerical data
    fscanf(infile, "%s", gene_name);

    // iterate over the columns of each row
    for (j = 0; j < cols; j++) {
      char element[50]; // used for reading in each element of the expression matrix
      if (fscanf(infile, "%s", element) != EOF) {
        // if this is a missing value and omission of missing values is enabled then
        // rewrite this value as MISSING_VALUE
        if (omit_na && strcmp(element, na_val) == 0) {
          data[i][j] = NAN;
        }
        else {
          if (do_log10) {
            data[i][j] = log10(atof(element));
          }
          else if (do_log2) {
            data[i][j] = log2(atof(element));
          }
          else if (do_log) {
            data[i][j] = log(atof(element));
          }
          else {
            data[i][j] = atof(element);
          }
        }
      }
      else {
        printf("Error: EOF reached early. Exiting.\n");
        exit(-1);
      }
    }
  }
  return data;
}

/**
 *
 */
int * init_histogram() {
  int * histogram = (int *) malloc(sizeof(int) * HIST_BINS + 1);

  int m; // integer used or looping

  // create and initialize the histogram array for the distribution of coefficients
  for (m = 0; m < HIST_BINS + 1; m++) {
    histogram[m] = 0;
  }
  return histogram;
}
/**
 *
 */
void print_histogram(int * histogram) {
  int m; // integer used or looping
  char outfilename[50];  // the output file name


  // output the correlation histogram
  sprintf(outfilename, "%s.mi.corrhist.txt", fileprefix);
  FILE * outfile = fopen(outfilename, "w");
  for (m = 0; m < HIST_BINS; m++) {
    fprintf(outfile, "%lf\t%d\n", 1.0 * m / (HIST_BINS), histogram[m]);
  }
  fclose(outfile);
}
/**
 *
 */
int calculate_MI(double ** data, int * histogram) {
  int i, j, k, m;   // integers for looping
  float one = 1.0;  // used for pairwise comparison of the same gene
  double mi;        // the mutual information value
  char outfilename[50];  // the output file name


  // output a maximum of ROWS_PER_OUTPUT_FILE rows and then start a new file
  int z = (rows - 1) / ROWS_PER_OUTPUT_FILE;
  int j_limit;

  // make sure the Pearson directory exists
  struct stat st = {0};
  if (stat("./MI", &st) == -1) {
      mkdir("./MI", 0700);
  }

  int total_comps = (rows * rows) / 2;
  int n_comps = 0;

  // each iteration of m is a new output file
  printf("Calculating mutual information...\n");
  for (m = 0; m <= z; m++) {

    // the output file will be located in the Pearson directory and named based on the input file info
    sprintf(outfilename, "./MI/%s.mi%d.bin", fileprefix, m);
    printf("Writing file %d of %d: %s... \n", m + 1, z + 1, outfilename);

    FILE * outfile = fopen(outfilename, "wb");

    // calculate the limit on the rows to output based on where we are in the calculation
    if (m < z) {
      j_limit = (m + 1) * ROWS_PER_OUTPUT_FILE;
    }
    else {
      j_limit = rows;
    }

    // output the size of the overall matrix
    fwrite(&rows, sizeof(rows), 1, outfile);

    // determine and output the number of rows that will be stored in the current file
    int i = j_limit - m * ROWS_PER_OUTPUT_FILE;
    fwrite(&i, sizeof(i), 1, outfile);

    for (j = m * ROWS_PER_OUTPUT_FILE; j < j_limit; j++) {
      // lower triangular symmetric matrix (stop when col# > row#)
      for (k = 0; k <= j; k++) {
        n_comps++;
        if (n_comps % 1000 == 0) {
          printf("Percent complete: %.2f%%\r", (n_comps/(float)total_comps)*100);
        }
        if (j == k) {
          // correlation of an element with itself is 1
          fwrite(&one, sizeof(one), 1, outfile);
        }
        else {
          // build the vectors for calculating MI
          double x[cols+1];
          double y[cols+1];
          int n = 0;
          for (i = 0; i < cols; i++) {
            // if either of these elements is missing then don't include the
            // elements from this sample
            if (data[j][i] == NAN || data[k][i] == NAN || isinf(data[j][i]) || isinf(data[k][i])) {
              continue;
            }
            x[n] = data[j][i];
            y[n] = data[k][i];
            n++;
          }

          // Calculate Mutual Information if we have enough observations.  We default the
          // MI value to NaN if we do not have the minimum number
          // of observations to do the calculation.  This is better than storing
          // a zero which indicates no correlation.  If we stored zero then
          // we'd be providing a false correlation as no correlation calculation
          // actually occurred.
          double mi = NAN;
          if (n >= min_obs) {
            mi = calculateMutualInformation(x, y, n);
          }
          fwrite(&mi, sizeof(float), 1, outfile);


          // if the historgram is turned on then store the value in the correct bin
          if (mi < 1 && mi > -1) {
            if (mi < 0) {
              mi = -mi;
            }
            histogram[(int)(mi * HIST_BINS)]++;
          }
        }
      }
    }
    fclose(outfile);
  }
  return 0;
}
/**
 * Main Subroutine
 */
int calculate_pearson(double ** data, int * histogram) {
  char gene_name[255];   // holds the gene or probe set name from the first column
  char outfilename[50];  // the output file name


  int i, j, k, m;   // integers for looping
  float one = 1.0;  // used for pairwise comparision of the same gene

  // output a maximum of ROWS_PER_OUTPUT_FILE rows and then start a new file
  int z = (rows - 1) / ROWS_PER_OUTPUT_FILE;
  int j_limit;

  // make sure the Pearson directory exists
  struct stat st = {0};
  if (stat("./Pearson", &st) == -1) {
      mkdir("./Pearson", 0700);
  }

  int total_comps = (rows * rows) / 2;
  int n_comps = 0;

  // each iteration of m is a new output file
  printf("Calculating correlations...\n");
  for (m = 0; m <= z; m++) {

    // the output file will be located in the Pearson directory and named based on the input file info
    sprintf(outfilename, "./Pearson/%s.pc%d.bin", fileprefix, m);
    printf("Writing file %d of %d: %s... \n", m + 1, z + 1, outfilename);

    FILE * outfile = fopen(outfilename, "wb");

    // calculate the limit on the rows to output based on where we are in the calculation
    if (m < z) {
      j_limit = (m + 1) * ROWS_PER_OUTPUT_FILE;
    }
    else {
      j_limit = rows;
    }

    // output the size of the overall matrix
    fwrite(&rows, sizeof(rows), 1, outfile);

    // determine and output the number of rows that will be stored in the current file
    i = j_limit - m * ROWS_PER_OUTPUT_FILE;
    fwrite(&i, sizeof(i), 1, outfile);

    for (j = m * ROWS_PER_OUTPUT_FILE; j < j_limit; j++) {
      // lower triangular symmetric matrix (stop when col# > row#) 
      for (k = 0; k <= j; k++) {
        n_comps++;
        if (n_comps % 1000 == 0) {
          printf("Percent complete: %.2f%%\r", (n_comps/(float)total_comps)*100);
        }
        if (j == k) {
          // correlation of an element with itself is 1
          fwrite(&one, sizeof(one), 1, outfile);
        }
        else {
          // build the vectors for calculating Pearson's
          double x[cols];
          double y[cols];
          int n = 0;
          for (i = 0; i < cols; i++) {
            // if either of these elements is missing then don't include the 
            // elements from this sample
            if (data[j][i] == NAN || data[k][i] == NAN || isinf(data[j][i]) || isinf(data[k][i])) {
              continue;
            }
            x[n] = data[j][i];
            y[n] = data[k][i];
            n++;
          }

          // Calculate Pearson's if we have enough observations.  We default the
          // Pearson value to NaN if we do not have the minimum number
          // of observations to do the calculation.  This is better than storing
          // a zero which indicates no correlation.  If we stored zero then
          // we'd be providing a false correlation as no correlation calculation
          // actually occurred.
          float pearson = NAN;
          if (n >= min_obs) {
            pearson = gsl_stats_correlation(x, 1, y, 1, n);
          }
          fwrite(&pearson, sizeof(float), 1, outfile);

          // if the historgram is turned on then store the value in the correct bin
          if (pearson < 1 && pearson > -1) {
            if (pearson < 0) {
              pearson = -pearson;
            }
            histogram[(int)(pearson * HIST_BINS)]++;
          }
        }
      }
    }
    fclose(outfile);
  }
  return 0;
}

/**
 *
 */
void print_usage() {
  printf("\n");
  printf("Usage: ./ccm [options]\n");
  printf("The list of required options:\n");
  printf("  --ematrix|-e The file name that contains the expression matrix.\n");
  printf("                 The rows must be genes or probe sets and columns are samples\n");
  printf("  --rows|-r    The number of lines in the input file including the header\n");
  printf("                 column if it exists\n");
  printf("  --cols|-c    The number of columns in the input file minus the first\n");
  printf("                 column that contains gene names\n");
  printf("  --method|-m  The correlation method to use. Supported methods include\n");
  printf("                 Pearson's correlation and Mutual Information. Provide\n");
  printf("                 either 'pc' or mi' as values respectively.\n");
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
  printf("  --hist        Provide this flag enable creation of correlation histogram.\n");
  printf("  --perf        Provide this flag to enable performance monitoring.\n");
  printf("  --headers     Provide this flag if the first line of the matrix contains\n");
  printf("                  headers.\n");
  printf("\n");
  printf("Note: Correlation values are set to NaN if there weren't enough observations\n");
  printf("to perform the calculation.\n");
}
