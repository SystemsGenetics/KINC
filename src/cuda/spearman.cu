
// #include "sort.cu"



/*!
 * Compute the Spearman correlation of a cluster in a pairwise data array.
 *
 * @param x
 * @param y
 * @param labels
 * @param sampleSize
 * @param cluster
 * @param minSamples
 * @param x_rank
 * @param y_rank
 */
__device__
float Spearman_computeCluster(
    const float *x,
    const float *y,
    const char *labels,
    int sampleSize,
    char cluster,
    int minSamples,
    float *x_rank,
    float *y_rank)
{
    // extract samples in pairwise cluster
    int n = 0;

    for ( int i = 0; i < sampleSize; ++i )
    {
        if ( labels[i] == cluster )
        {
            x_rank[n] = x[i];
            y_rank[n] = y[i];
            ++n;
        }
    }

    // get power of 2 size
    int N_pow2 = nextPower2(sampleSize);

    for ( int i = n; i < N_pow2; ++i )
    {
        x_rank[i] = INFINITY;
        y_rank[i] = INFINITY;
    }

    // compute correlation only if there are enough samples
    float result = NAN;

    if ( n >= minSamples )
    {
        // compute rank of x
        bitonicSortFF(N_pow2, x_rank, y_rank);
        computeRank(x_rank, n);

        // compute rank of y
        bitonicSortFF(N_pow2, y_rank, x_rank);
        computeRank(y_rank, n);

        // compute correlation of rank arrays
        float sumx = 0;
        float sumy = 0;
        float sumx2 = 0;
        float sumy2 = 0;
        float sumxy = 0;

        for ( int i = 0; i < n; ++i )
        {
            float x_i = x_rank[i];
            float y_i = y_rank[i];

            sumx += x_i;
            sumy += y_i;
            sumx2 += x_i * x_i;
            sumy2 += y_i * y_i;
            sumxy += x_i * y_i;
        }

        result = (n*sumxy - sumx*sumy) / sqrt((n*sumx2 - sumx*sumx) * (n*sumy2 - sumy*sumy));
    }

    return result;
}



/*!
 * Compute the correlation of each cluster in a pairwise data array. The data array
 * should only contain the clean samples that were extracted from the expression
 * matrix, while the labels should contain all samples.
 *
 * @param x
 * @param y
 * @param sampleSize
 * @param clusterSize
 * @param labels
 * @param minSamples
 * @param x_rank
 * @param y_rank
 * @param correlations
 */
__device__
void Spearman_compute(
    const float *x,
    const float *y,
    int sampleSize,
    char clusterSize,
    const char *labels,
    int minSamples,
    float *x_rank,
    float *y_rank,
    float *correlations)
{
    for ( char k = 0; k < clusterSize; ++k )
    {
        correlations[k] = Spearman_computeCluster(
            x, y,
            labels,
            sampleSize,
            k,
            minSamples,
            x_rank,
            y_rank
        );
    }
}
