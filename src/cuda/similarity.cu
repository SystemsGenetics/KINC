
// #include "fetchpair.cu"
// #include "outlier.cu"
// #include "gmm.cu"
// #include "pearson.cu"
// #include "spearman.cu"



enum ClusteringMethod
{
    ClusteringMethod_None
    ,ClusteringMethod_GMM
};



enum CorrelationMethod
{
    CorrelationMethod_Pearson
    ,CorrelationMethod_Spearman
};



/*!
 * Compute the similarity measure for the given block of pairwise indices.
 *
 * @param clusMethod
 * @param corrMethod
 * @param removePreOutliers
 * @param removePostOutliers
 * @param numPairs
 * @param expressions
 * @param sampleSize
 * @param in_index
 * @param minExpression
 * @param minSamples
 * @param minClusters
 * @param maxClusters
 * @param criterion
 * @param work_x
 * @param work_y
 * @param work_gmm_data
 * @param work_gmm_labels
 * @param work_gmm_pi
 * @param work_gmm_mu
 * @param work_gmm_sigma
 * @param work_gmm_sigmaInv
 * @param work_gmm_normalizer
 * @param work_gmm_MP
 * @param work_gmm_counts
 * @param work_gmm_logpi
 * @param work_gmm_gamma
 * @param out_K
 * @param out_labels
 * @param out_correlations
 */
__global__
void Similarity_compute(
    ClusteringMethod  clusMethod,
    CorrelationMethod corrMethod,
    bool              removePreOutliers,
    bool              removePostOutliers,
    int               numPairs,
    const float *     expressions,
    int               sampleSize,
    const int2 *      in_index,
    int               minExpression,
    int               minSamples,
    char              minClusters,
    char              maxClusters,
    Criterion         criterion,
    float *           work_x,
    float *           work_y,
    Vector2 *         work_gmm_data,
    char *            work_gmm_labels,
    float *           work_gmm_pi,
    Vector2 *         work_gmm_mu,
    Matrix2x2 *       work_gmm_sigma,
    Matrix2x2 *       work_gmm_sigmaInv,
    float *           work_gmm_normalizer,
    Vector2 *         work_gmm_MP,
    int *             work_gmm_counts,
    float *           work_gmm_logpi,
    float *           work_gmm_gamma,
    char *            out_K,
    char *            out_labels,
    float *           out_correlations)
{
    int offset = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;

    for ( int i = offset; i < numPairs; i += stride )
    {
        // initialize variables
        int N_pow2 = nextPower2(sampleSize);

        int2 index = in_index[i];
        const float *x = &expressions[index.x * sampleSize];
        const float *y = &expressions[index.y * sampleSize];
        float *x_sorted = &work_x[i * N_pow2];
        float *y_sorted = &work_y[i * N_pow2];

        GMM gmm = {
            &work_gmm_data[i * sampleSize],
            &work_gmm_labels[i * sampleSize],
            &work_gmm_pi[i * maxClusters],
            &work_gmm_mu[i * maxClusters],
            &work_gmm_sigma[i * maxClusters],
            &work_gmm_sigmaInv[i * maxClusters],
            &work_gmm_normalizer[i * maxClusters],
            &work_gmm_MP[i * maxClusters],
            &work_gmm_counts[i * maxClusters],
            &work_gmm_logpi[i * maxClusters],
            &work_gmm_gamma[i * maxClusters * sampleSize],
            0,
            0
        };

        char *p_K = &out_K[i];
        char *labels = &out_labels[i * sampleSize];
        float *correlations = &out_correlations[i * maxClusters];

        // fetch pairwise input data
        int numSamples = fetchPair(
            x, y,
            sampleSize,
            minExpression,
            labels
        );

        // remove pre-clustering outliers
        if ( removePreOutliers )
        {
            numSamples = removeOutliers(
                x, y,
                sampleSize,
                numSamples,
                labels,
                1,
                -7,
                x_sorted,
                y_sorted
            );
        }

        // compute clusters
        char K = 1;

        if ( clusMethod == ClusteringMethod_GMM )
        {
            K = GMM_compute(
                &gmm,
                x, y,
                sampleSize,
                numSamples,
                labels,
                minSamples,
                minClusters,
                maxClusters,
                criterion
            );
        }

        // remove post-clustering outliers
        if ( removePostOutliers )
        {
            numSamples = removeOutliers(
                x, y,
                sampleSize,
                numSamples,
                labels,
                K,
                -8,
                x_sorted,
                y_sorted
            );
        }

        // compute correlations
        if ( corrMethod == CorrelationMethod_Pearson )
        {
            Pearson_compute(
                x, y,
                sampleSize,
                K,
                labels,
                minSamples,
                correlations
            );
        }
        else if ( corrMethod == CorrelationMethod_Spearman )
        {
            Spearman_compute(
                x, y,
                sampleSize,
                K,
                labels,
                minSamples,
                x_sorted,
                y_sorted,
                correlations
            );
        }

        // save number of clusters
        *p_K = K;
    }
}
