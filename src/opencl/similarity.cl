
// #include "fetchpair.cl"
// #include "outlier.cl"
// #include "gmm.cl"
// #include "pearson.cl"
// #include "spearman.cl"



typedef enum
{
    ClusteringMethod_None
    ,ClusteringMethod_GMM
} ClusteringMethod;



typedef enum
{
    CorrelationMethod_Pearson
    ,CorrelationMethod_Spearman
} CorrelationMethod;



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
 * @param maxExpression
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
__kernel
void Similarity_compute(
    ClusteringMethod       clusMethod,
    CorrelationMethod      corrMethod,
    int                    removePreOutliers,
    int                    removePostOutliers,
    int                    numPairs,
    __global const float * expressions,
    int                    sampleSize,
    __global const int2 *  in_index,
    float                  minExpression,
    float                  maxExpression,
    int                    minSamples,
    char                   minClusters,
    char                   maxClusters,
    Criterion              criterion,
    __global float *       work_x,
    __global float *       work_y,
    __global Vector2 *     work_gmm_data,
    __global char *        work_gmm_labels,
    __global float *       work_gmm_pi,
    __global Vector2 *     work_gmm_mu,
    __global Matrix2x2 *   work_gmm_sigma,
    __global Matrix2x2 *   work_gmm_sigmaInv,
    __global float *       work_gmm_normalizer,
    __global Vector2 *     work_gmm_MP,
    __global int *         work_gmm_counts,
    __global float *       work_gmm_logpi,
    __global float *       work_gmm_gamma,
    __global char *        out_K,
    __global char *        out_labels,
    __global float *       out_correlations)
{
    int i = get_global_id(0);

    if ( i >= numPairs )
    {
        return;
    }

    // initialize variables
    int N_pow2 = nextPower2(sampleSize);

    int2 index = in_index[i];
    __global const float *x = &expressions[index.x * sampleSize];
    __global const float *y = &expressions[index.y * sampleSize];
    __global float *x_sorted = &work_x[i * N_pow2];
    __global float *y_sorted = &work_y[i * N_pow2];

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

    __global char *p_K = &out_K[i];
    __global char *labels = &out_labels[i * sampleSize];
    __global float *correlations = &out_correlations[i * maxClusters];

    // fetch pairwise input data
    int numSamples = fetchPair(
        x, y,
        sampleSize,
        minExpression,
        maxExpression,
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
