#include "PearsonSimilarity.h"



/**
 *
 */
PearsonSimilarity::PearsonSimilarity(PairWiseSet * pws, int min_obs)
  :PairWiseSimilarity(pws, min_obs) {

  strcpy(this->type, "pc");
}
/**
 *
 */
PearsonSimilarity::PearsonSimilarity(PairWiseSet * pws, int min_obs, int * samples)
  :PairWiseSimilarity(pws, min_obs, samples) {

  strcpy(this->type, "pc");
}

PearsonSimilarity::~PearsonSimilarity() {

}

/**
 * Performs Pearson correlation on two arrays.
 *
 * @param float *a
 * @param float *b
 * @param int n
 */
void PearsonSimilarity::run() {
  // Make sure we have the correct number of observations before performing
  // the comparision.
  if (this->n >= this->min_obs) {
    score = gsl_stats_float_correlation(this->a, 1, this->b, 1, this->n);
  }
  else {
    score = NAN;
  }
}

/**
 * Calculates the Pearson correlation matrix
 *
 * @param CCMParameters params
 *   An instance of the CCM parameters struct
 * @param float ** data
 *   A pointer to a two dimensional array of doubles containing the
 *   expression values
 * @param int * histogram
 *   An integer array created by the init_histogram function and
 *   used for building the histogram.
 *
 */
/*void calculatePearson(CCMParameters params, float ** data, int * histogram) {
  char outfilename[1024];  // the output file name

  int i, j, k;      // integers for looping
  float one = 1.0;  // used for pairwise comparision of the same gene
  int num_bins;     // the number of binary files needed to store the matrix
  int curr_bin;     // the current binary file number
  int bin_rows;     // holds the number of rows in the file
  int total_comps;  // the total number of pair-wise comparisions to be made
  int n_comps;      // the number of comparisions completed during looping

  // calculate the number of binary files needed to store the similarity matrix
  num_bins = (params.rows - 1) / ROWS_PER_OUTPUT_FILE;

  // make sure the Pearson directory exists
  struct stat st = {0};
  if (stat("./Pearson", &st) == -1) {
      mkdir("./Pearson", 0700);
  }

  total_comps = (params.rows * params.rows) / 2;
  n_comps = 0;

  // each iteration of m is a new output file
  printf("Calculating correlations...\n");
  for (curr_bin = 0; curr_bin <= num_bins; curr_bin++) {

    // calculate the limit on the rows to output based on where we are in the calculation
    if (curr_bin < num_bins) {
      bin_rows = (curr_bin + 1) * ROWS_PER_OUTPUT_FILE;
    }
    else {
      bin_rows = params.rows;
    }
    int num_cols = bin_rows - (curr_bin * ROWS_PER_OUTPUT_FILE);

    // the output file will be located in the Pearson directory and named based on the input file info
    sprintf(outfilename, "./Pearson/%s.pc%d.bin", params.fileprefix, curr_bin);
    printf("Writing file %d of %d: %s... \n", curr_bin + 1, num_bins + 1, outfilename);
    FILE * outfile = fopen(outfilename, "wb");

    // write the size of the matrix
    fwrite(&params.rows, sizeof(params.rows), 1, outfile);
    // write the number of lines in this file
    fwrite(&num_cols, sizeof(num_cols), 1, outfile);

    // iterate through the genes that belong in this file
    for (j = curr_bin * ROWS_PER_OUTPUT_FILE; j < bin_rows; j++) {
      // iterate through all the genes up to j (only need lower triangle)
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
          float x[params.cols];
          float y[params.cols];
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
            pearson = gsl_stats_float_correlation(x, 1, y, 1, n);
            //printf("(%d, %d) = %f (%d values).\n", j, k, pearson, n);
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
}*/
