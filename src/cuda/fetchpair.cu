
// #include "linalg.cu"



/*!
 * Compute the initial labels for a gene pair in an expression matrix. Samples
 * with missing values and samples that are outside the expression thresholds are
 * labeled as such, all other samples are labeled as cluster 0. The number of
 * clean samples is returned.
 *
 * @param x
 * @param y
 * @param sampleSize
 * @param minExpression
 * @param maxExpression
 * @param labels
 */
__device__
int fetchPair(
    const float *x,
    const float *y,
    int sampleSize,
    float minExpression,
    float maxExpression,
    char *labels)
{
    // label the pairwise samples
    int N = 0;

    for ( int i = 0; i < sampleSize; ++i )
    {
        // label samples with missing values
        if ( isnan(x[i]) || isnan(y[i]) )
        {
            labels[i] = -9;
        }

        // label samples which are below the minimum expression threshold
        else if ( x[i] < minExpression || y[i] < minExpression )
        {
            labels[i] = -6;
        }

        // label samples which are above the maximum expression threshold
        else if ( x[i] > maxExpression || y[i] > maxExpression )
        {
            labels[i] = -6;
        }

        // label any remaining samples as cluster 0
        else
        {
            N++;
            labels[i] = 0;
        }
    }

    // return number of clean samples
    return N;
}
