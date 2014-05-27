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

#include <gsl/gsl_statistics.h>

/**
 * Constants
 */

// a global variable for the number of rows in each output file
int ROWS_PER_OUTPUT_FILE = 10000;

// the number of bins in the historgram of correlation values
int HIST_BINS = 100;

// specify the value to be used as a missing value
int MISSING_VALUE = -1;


/**
 * Main Subroutine
 */
int main(int argc, char *argv[]) {

  // variables used for file I/O
  FILE *infile, *outfile;// pointers to the input and output files
  char gene_name[255];   // holds the gene or probeset name from the first column
  char outfilename[50];  // the output file name
  char *infilename;      // the input file name
  char fileprefix[1024]; // the input filename without the prefix 

  // variables specified by the user on the command-line
  int rows;              // the number of rows in the expression matrix
  int cols;              // the number of cols in the expression matrix
  int perf = 0;          // indicates if performance monitoring should be enabled
  int omit_na = 0;       // indicates if missing values should be ignored
  int min_obs = 30;      // the minimum number of observations to calculate correlation
  int do_hist = 0;       // set to 1 to enable creation of correlation histogram
  char *na_val;          // specifies the value that represents a missing value
  char func[10];         // specifies the transformation function: log2, none
  int do_log2;           // set to 1 to perform log2 transformation
  int headers = 0;       // set to 1 if the first line contains headers

  // variables used for calculating Pearson's correlation
  int i, j, k, m;
  double **data;          // holds the expressionn matrix
  float one = 1.0;       // used for pairwise comparision of the same gene

  // variables used for timing of the software
  time_t start_time, end_time;


  // make sure we have between 4 and 5 incoming arguments
  if(argc < 4 || argc > 11) {
    printf("Usage: ./ccm <ematrix> <rows> <cols> [<omit_na> <na_val> <hist> <perf> <headers>]\n");
    printf("  <ematrix>: the file name that contains the expression matrix. The rows must be genes or probesets and columns are samples\n");
    printf("  <rows>:    the number of lines in the input file minus the header column if it exists\n");
    printf("  <cols>:    the number of columns in the input file minus the first column that contains gene names\n");
    printf("  <omit_na>: set to 1 to ignore missing values. Defaults to 0.\n");
    printf("  <na_val>:  a string representing the missing values in the input file (e.g. NA or 0.000)\n");
    printf("  <min_obs>: the minimum number of observations (after missing values removed) that must be present to perform correlation. Default is 30\n");
    printf("  <func>:    a transformation function to apply to elements of the ematrix. Values include: log2 or none. Default is none\n");
    printf("  <hist>:    set to 1 to enable creation of correlation historgram. Defaults to 0.\n");
    printf("  <perf>:    set to 1 to enable performance monitoring. Defaults to 0\n");
    printf("  <headers>: set to 1 if the first line contains headers. Defaults to 0\n");
    printf("Note: Correlation value is set to NaN if there weren't enough observations to perform the calculation.\n");
    exit(-1);
  }

  // get the incoming arguments
  infilename = argv[1];
  rows = atoi(argv[2]);  // the number of rows in the expression matrix
  cols = atoi(argv[3]);  // the number of cols in the expression matrix

  // get the incoming arguments related to missing values
  if (argc >= 5) {
    omit_na = atoi(argv[4]);
  }
  if (argc >= 6) {
    na_val = argv[5];
  }
  if (argc >= 7) {
    min_obs = atoi(argv[6]);
  }
  strcpy(func, "none");
  do_log2 = 0;
  if (argc >= 8 && strcmp(argv[7], "log2") == 0) {
    strcpy(func, "log2");
    do_log2 = 1;
  }
  if (argc >= 9 && atoi(argv[8]) == 1) {
    do_hist = 1;
  }

  // enable performance monitoring if 10th argument is set
  if (argc >= 10 && atoi(argv[9]) == 1) {
    perf = 1;
  }

  // skip header line if 11th argument is set
  if (argc == 11 && atoi(argv[10]) == 1) {
    headers = 1;
  }

  // if performance monitoring is enabled the set the start time
  if (perf) {
    time(&start_time);
  }

  // remove the path and extension from the filename
  strcpy(fileprefix, basename(infilename));
  char *p = rindex(fileprefix, '.');
  p[0] = 0;

  // allocate the data array for storing the input expression matrix
  data = (double**) malloc(sizeof(double *) * rows);
  for (i = 0; i < rows; i++) {
    data[i] = malloc(sizeof(double) * cols);
  }

  // read in expression matrix from the user-specified file
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

    for (j = 0; j < cols; j++) {
      char element[50]; // used for reading in each element of the expression matrix
      if (fscanf(infile, "%s", element) != EOF) {
        // if this is a missing value and omission of missing values is enabled then 
        // rewrite this value as MISSING_VALUE
        if (omit_na && strcmp(element, na_val) == 0) {
          data[i][j] = MISSING_VALUE;
        }
        else {
          if (do_log2) {
            data[i][j] = log2(atof(element));
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

  // output a maximum of ROWS_PER_OUTPUT_FILE rows and then start a new file
  int z = (rows - 1) / ROWS_PER_OUTPUT_FILE;
  int j_limit;

  // create and initialize the histogram array for the distribution of coefficients
  int histogram[HIST_BINS + 1];
  if (do_hist) {
    for (m = 0; m < HIST_BINS + 1; m++) {
      histogram[m] = 0;
    }
  }

  // make sure the Pearson directory exists
  struct stat st = {0};
  if (stat("./Pearson", &st) == -1) {
      mkdir("./Pearson", 0700);
  }

  // each iteration of m is a new output file
  printf("Calculating correlations...\n");
  for (m = 0; m <= z; m++) {

    // the output file will be located in the Pearson directory and named based on the input file info
    sprintf(outfilename, "./Pearson/%s.pc%d.bin", fileprefix, m);
    printf("Writing file %d of %d: %s... \n", m + 1, z + 1, outfilename);

    outfile = fopen(outfilename, "wb");

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
            if (data[j][i] == MISSING_VALUE || data[k][i] == MISSING_VALUE) {
              continue;
            }
            x[n] = data[j][i];
            y[n] = data[k][i];
            n++;
          }

          // calculate Pearsons if we have enough observations.  We default the
          // pearson value to NaN if we do not have the minimum number
          // of observations to do the calculation.  This is better than storing
          // a zero which indicates no correlation.  If we stored zero then
          // we'd be providing a false correlation as no correlation calculation
          // actually occured.
          float pearson = NAN;
          if (n >= min_obs) {
            pearson = gsl_stats_correlation(x, 1, y, 1, n);
          }
          fwrite(&pearson, sizeof(float), 1, outfile);

          // if the historgram is turned on then store the value in the correct bin
          if (do_hist && pearson < 1 && pearson > -1) {
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

  // if do_hist is not zero then output the correlation histogram
  if (do_hist) {
    outfile = fopen("corrhist.txt", "w");
    for (m = 0; m < HIST_BINS; m++) {
      fprintf(outfile, "%lf\t%d\n", 1.0 * m / (HIST_BINS), histogram[m]);
    }
    fclose(outfile);
  }

  // if performance monitoring is enabled then write the timing data
  if(perf) {
    time(&end_time);
    outfile = fopen("timingdata.txt", "a");
    fprintf(outfile, "CCM Runtime with %d x %d %s input dataset: %.2lf min\n", rows, cols, infilename, difftime(end_time, start_time)/60.0);
  }

  printf("Done.\n");
  return 0;
}
