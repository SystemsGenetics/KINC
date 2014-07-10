#include "pearson.h"

/**
 * Calculates the Pearson correlation matrix
 *
 * @param CCMParameters params
 *   An instance of the CCM parameters struct
 * @param double ** data
 *   A pointer to a two dimensional array of doubles containing the
 *   expression values
 * @param int * histogram
 *   An integer array created by the init_histogram function and
 *   used for building the histogram.
 *
 */
void calculate_pearson(CCMParameters params, double ** data, int * histogram) {
  char outfilename[50];  // the output file name


  int i, j, k, m;   // integers for looping
  float one = 1.0;  // used for pairwise comparision of the same gene

  // output a maximum of ROWS_PER_OUTPUT_FILE rows and then start a new file
  int z = (params.rows - 1) / ROWS_PER_OUTPUT_FILE;
  int j_limit;

  // make sure the Pearson directory exists
  struct stat st = {0};
  if (stat("./Pearson", &st) == -1) {
      mkdir("./Pearson", 0700);
  }

  int total_comps = (params.rows * params.rows) / 2;
  int n_comps = 0;

  // each iteration of m is a new output file
  printf("Calculating correlations...\n");
  for (m = 0; m <= z; m++) {

    // the output file will be located in the Pearson directory and named based on the input file info
    sprintf(outfilename, "./Pearson/%s.pc%d.bin", params.fileprefix, m);
    printf("Writing file %d of %d: %s... \n", m + 1, z + 1, outfilename);

    FILE * outfile = fopen(outfilename, "wb");

    // calculate the limit on the rows to output based on where we are in the calculation
    if (m < z) {
      j_limit = (m + 1) * ROWS_PER_OUTPUT_FILE;
    }
    else {
      j_limit = params.rows;
    }

    // output the size of the overall matrix
    fwrite(&params.rows, sizeof(params.rows), 1, outfile);

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
          double x[params.cols];
          double y[params.cols];
          int n = 0;
          for (i = 0; i < params.cols; i++) {
            // if either of these elements is missing then don't include the
            // elements from this sample
            if (isnan(data[j][i]) || isnan(data[k][i]) || isinf(data[j][i]) || isinf(data[k][i])) {
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
          if (n >= params.min_obs) {
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
}
