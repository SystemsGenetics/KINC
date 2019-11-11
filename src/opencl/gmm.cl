
// #include "fetchpair.cl"
// #include "linalg.cl"
// #include "outlier.cl"



typedef struct
{
    __global Vector2 *   data;
    __global char *      labels;
    __global float *     pi;
    __global Vector2 *   mu;
    __global Matrix2x2 * sigma;
    __global Matrix2x2 * sigmaInv;
    __global float *     normalizer;
    __global Vector2 *   MP;
    __global int *       counts;
    __global float *     logpi;
    __global float *     gamma;
    float                logL;
    float                entropy;
} GMM;



/*!
 * Implementation of rand(), taken from POSIX example.
 *
 * @param state
 */
int myrand(ulong *state)
{
    *state = (*state) * 1103515245 + 12345;
    return ((unsigned)((*state)/65536) % 32768);
}



/*!
 * Initialize mixture components.
 *
 * @param gmm
 * @param X
 * @param N
 * @param K
 */
void GMM_initializeComponents(
    GMM *gmm,
    __global const Vector2 *X,
    int N,
    int K)
{
    // initialize random state
    ulong state = 1;

    // initialize each mixture component
    for ( int k = 0; k < K; ++k )
    {
        // initialize mixture weight to uniform distribution
        gmm->pi[k] = 1.0f / K;

        // initialize mean to a random sample from X
        int i = myrand(&state) % N;

        gmm->mu[k] = X[i];

        // initialize covariance to identity matrix
        matrixInitIdentity(&gmm->sigma[k]);
    }
}



/*!
 * Pre-compute the precision matrix and normalizer term for each
 * mixture component.
 *
 * @param gmm
 * @param K
 */
bool GMM_prepareComponents(GMM *gmm, int K)
{
    const int D = 2;
    bool success = true;

    for ( int k = 0; k < K; ++k )
    {
        // compute precision (inverse of covariance)
        float det;
        matrixInverse(&gmm->sigma[k], &gmm->sigmaInv[k], &det);

        // compute normalizer term for multivariate normal distribution
        gmm->normalizer[k] = -0.5f * (D * log(2.0f * M_PI) + log(det));

        // return failure if matrix inverse failed
        success = success && (det > 0) && !isnan(det);
    }

    return success;
}



/*!
 * Initialize the mean of each component in the mixture model using k-means
 * clustering.
 *
 * @param gmm
 * @param X
 * @param N
 * @param K
 */
void GMM_initializeMeans(GMM *gmm, __global const Vector2 *X, int N, int K)
{
    const int MAX_ITERATIONS = 20;
    const float TOLERANCE = 1e-3f;
    float diff = 0;

    // initialize workspace
    __global Vector2 *MP = gmm->MP;
    __global int *counts = gmm->counts;

    for ( int t = 0; t < MAX_ITERATIONS && diff > TOLERANCE; ++t )
    {
        // compute mean and sample count for each component
        for ( int k = 0; k < K; ++k )
        {
            vectorInitZero(&MP[k]);
            counts[k] = 0;
        }

        for ( int i = 0; i < N; ++i )
        {
            // determine the component mean which is nearest to x_i
            float min_dist = INFINITY;
            int min_k = 0;
            for ( int k = 0; k < K; ++k )
            {
                float dist = vectorDiffNorm(&X[i], &gmm->mu[k]);
                if ( min_dist > dist )
                {
                    min_dist = dist;
                    min_k = k;
                }
            }

            // update mean and sample count
            vectorAdd(&MP[min_k], &X[i]);
            ++counts[min_k];
        }

        // scale each mean by its sample count
        for ( int k = 0; k < K; ++k )
        {
            vectorScale(&MP[k], 1.0f / counts[k]);
        }

        // compute the total change of all means
        diff = 0;
        for ( int k = 0; k < K; ++k )
        {
            diff += vectorDiffNorm(&MP[k], &gmm->mu[k]);
        }
        diff /= K;

        // update component means
        for ( int k = 0; k < K; ++k )
        {
            gmm->mu[k] = MP[k];
        }
    }
}



/*!
 * Perform the expectation step of the EM algorithm. In this step we compute
 * gamma, the posterior probabilities for each component in the mixture model
 * and each sample in X, as well as the log-likelihood of the model.
 *
 * First we compute the probability density function of the multivariate normal
 * distribution conditioned on a single component for each point in X:
 *
 *   P(x|k) = exp(-0.5 * (x - mu_k)^T Sigma_k^-1 (x - mu_k)) / sqrt((2pi)^d det(Sigma_k))
 *
 * We actually use the log-probability to avoid numerical instability:
 *
 *   log(P(x|k)) = -0.5 * (x - mu_k)^T Sigma_k^-1 (x - mu_k) - 0.5 * (d * log(2pi) + log(det(Sigma_k)))
 *
 * Then we can compute gamma, the posterior probability matrix, and log(L), the
 * log-likelihood:
 *
 *   log(p(x_i)) = a + log(sum(exp(log(pi_k) + log(P(x_i|k))) - a))
 *
 *   gamma_ki = exp(log(pi_k) + log(P(x_i|k)) - log(p(x_i)))
 *
 *   log(L) = sum(log(p(x_i)))
 *
 * @param gmm
 * @param X
 * @param N
 * @param K
 */
float GMM_computeEStep(GMM *gmm, __global const Vector2 *X, int N, int K)
{
    // compute logpi
    for ( int k = 0; k < K; ++k )
    {
        gmm->logpi[k] = log(gmm->pi[k]);
    }

    // compute the log-probability for each component and each point in X
    __global float *logProb = gmm->gamma;

    for ( int k = 0; k < K; ++k )
    {
        for ( int i = 0; i < N; ++i )
        {
            // compute xm = (x - mu)
            Vector2 xm = X[i];
            vectorSubtract(&xm, &gmm->mu[k]);

            // compute Sxm = Sigma^-1 xm
            Vector2 Sxm;
            matrixProduct(&gmm->sigmaInv[k], &xm, &Sxm);

            // compute xmSxm = xm^T Sigma^-1 xm
            float xmSxm = vectorDot(&xm, &Sxm);

            // compute log(P) = normalizer - 0.5 * xm^T * Sigma^-1 * xm
            logProb[k * N + i] = gmm->normalizer[k] - 0.5f * xmSxm;
        }
    }

    // compute gamma and log-likelihood
    float logL = 0;

    for ( int i = 0; i < N; ++i )
    {
        // compute a = argmax(logpi_k + logProb_ki, k)
        float maxArg = -INFINITY;
        for ( int k = 0; k < K; ++k )
        {
            float arg = gmm->logpi[k] + logProb[k * N + i];
            if ( maxArg < arg )
            {
                maxArg = arg;
            }
        }

        // compute logpx
        float sum = 0;
        for ( int k = 0; k < K; ++k )
        {
            sum += exp(gmm->logpi[k] + logProb[k * N + i] - maxArg);
        }

        float logpx = maxArg + log(sum);

        // compute gamma_ki
        for ( int k = 0; k < K; ++k )
        {
            gmm->gamma[k * N + i] += gmm->logpi[k] - logpx;
            gmm->gamma[k * N + i] = exp(gmm->gamma[k * N + i]);
        }

        // update log-likelihood
        logL += logpx;
    }

    // return log-likelihood
    return logL;
}



/*!
 * Perform the maximization step of the EM algorithm. In this step we update the
 * parameters of the the mixture model using gamma, which is computed during the
 * expectation step:
 *
 *   n_k = sum(gamma_ki)
 *
 *   pi_k = n_k / N
 *
 *   mu_k = sum(gamma_ki * x_i)) / n_k
 *
 *   Sigma_k = sum(gamma_ki * (x_i - mu_k) * (x_i - mu_k)^T) / n_k
 *
 * @param gmm
 * @param X
 * @param N
 * @param K
 */
void GMM_computeMStep(GMM *gmm, __global const Vector2 *X, int N, int K)
{
    for ( int k = 0; k < K; ++k )
    {
        // compute n_k = sum(gamma_ki)
        float n_k = 0;

        for ( int i = 0; i < N; ++i )
        {
            n_k += gmm->gamma[k * N + i];
        }

        // update mixture weight
        gmm->pi[k] = n_k / N;

        // update mean
        Vector2 mu;

        vectorInitZero(&mu);

        for ( int i = 0; i < N; ++i )
        {
            vectorAddScaled(&mu, gmm->gamma[k * N + i], &X[i]);
        }

        vectorScale(&mu, 1.0f / n_k);

        gmm->mu[k] = mu;

        // update covariance matrix
        Matrix2x2 sigma;

        matrixInitZero(&sigma);

        for ( int i = 0; i < N; ++i )
        {
            // compute xm = (x_i - mu_k)
            Vector2 xm = X[i];
            vectorSubtract(&xm, &mu);

            // compute Sigma_ki = gamma_ki * (x_i - mu_k) (x_i - mu_k)^T
            matrixAddOuterProduct(&sigma, gmm->gamma[k * N + i], &xm);
        }

        matrixScale(&sigma, 1.0f / n_k);

        gmm->sigma[k] = sigma;
    }
}



/*!
 * Compute the cluster labels of a dataset using gamma:
 *
 *   y_i = argmax(gamma_ki, k)
 *
 * @param gamma
 * @param N
 * @param K
 * @param labels
 */
void GMM_computeLabels(
    __global const float *gamma,
    int N,
    int K,
    __global char *labels)
{
    for ( int i = 0; i < N; ++i )
    {
        // determine the value k for which gamma_ki is highest
        int max_k = -1;
        float max_gamma = -INFINITY;

        for ( int k = 0; k < K; ++k )
        {
            if ( max_gamma < gamma[k * N + i] )
            {
                max_k = k;
                max_gamma = gamma[k * N + i];
            }
        }

        // assign x_i to cluster k
        labels[i] = max_k;
    }
}



/*!
 * Compute the entropy of the mixture model for a dataset using gamma
 * and the given cluster labels:
 *
 *   E = -sum(sum(z_ki * log(gamma_ki))), z_ki = (y_i == k)
 *
 * @param gamma
 * @param N
 * @param labels
 */
float GMM_computeEntropy(
    __global const float *gamma,
    int N,
    __global const char *labels)
{
    float E = 0;

    for ( int i = 0; i < N; ++i )
    {
        int k = labels[i];

        E -= log(gamma[k * N + i]);
    }

    return E;
}



/*!
 * Fit the mixture model to a pairwise data array and compute the output cluster
 * labels for the data. The data array should only contain clean samples.
 *
 * @param gmm
 * @param X
 * @param N
 * @param K
 * @param labels
 */
bool GMM_fit(
    GMM *gmm,
    __global const Vector2 *X,
    int N,
    int K,
    __global char *labels)
{
    // initialize mixture components
    GMM_initializeComponents(gmm, X, N, K);

    // initialize means with k-means
    GMM_initializeMeans(gmm, X, N, K);

    // run EM algorithm
    const int MAX_ITERATIONS = 100;
    const float TOLERANCE = 1e-8f;
    float prevLogL = -INFINITY;
    float currLogL = -INFINITY;

    for ( int t = 0; t < MAX_ITERATIONS; ++t )
    {
        // pre-compute precision matrix and normalizer term for each mixture component
        bool success = GMM_prepareComponents(gmm, K);

        // return failure if matrix inverse failed
        if ( !success )
        {
            return false;
        }

        // perform E step
        prevLogL = currLogL;
        currLogL = GMM_computeEStep(gmm, X, N, K);

        // check for convergence
        if ( fabs(currLogL - prevLogL) < TOLERANCE )
        {
            break;
        }

        // perform M step
        GMM_computeMStep(gmm, X, N, K);
    }

    // save outputs
    gmm->logL = currLogL;
    GMM_computeLabels(gmm->gamma, N, K, labels);
    gmm->entropy = GMM_computeEntropy(gmm->gamma, N, labels);

    return true;
}



typedef enum
{
    AIC,
    BIC,
    ICL
} Criterion;



/*!
 * Compute the Akaike Information Criterion of a Gaussian mixture model.
 *
 * @param K
 * @param D
 * @param logL
 */
float GMM_computeAIC(int K, int D, float logL)
{
    int p = K * (1 + D + D * D);

    return 2 * p - 2 * logL;
}



/*!
 * Compute the Bayesian Information Criterion of a Gaussian mixture model.
 *
 * @param K
 * @param D
 * @param logL
 * @param N
 */
float GMM_computeBIC(int K, int D, float logL, int N)
{
    int p = K * (1 + D + D * D);

    return log((float) N) * p - 2 * logL;
}



/*!
 * Compute the Integrated Completed Likelihood of a Gaussian mixture model.
 *
 * @param K
 * @param D
 * @param logL
 * @param N
 * @param E
 */
float GMM_computeICL(int K, int D, float logL, int N, float E)
{
    int p = K * (1 + D + D * D);

    return log((float) N) * p - 2 * logL + 2 * E;
}



/*!
 * Determine the number of clusters in a pairwise data array. Several sub-models,
 * each one having a different number of clusters, are fit to the data and the
 * sub-model with the best criterion value is selected.
 *
 * @param gmm
 * @param x
 * @param y
 * @param sampleSize
 * @param numSamples
 * @param labels
 * @param minSamples
 * @param minClusters
 * @param maxClusters
 * @param criterion
 */
char GMM_compute(
    GMM *gmm,
    __global const float *x,
    __global const float *y,
    int sampleSize,
    int numSamples,
    __global char *labels,
    int minSamples,
    char minClusters,
    char maxClusters,
    Criterion criterion)
{
    // perform clustering only if there are enough samples
    char bestK = 0;

    if ( numSamples >= minSamples )
    {
        // extract clean samples from data array
        for ( int i = 0, j = 0; i < sampleSize; ++i )
        {
            if ( labels[i] >= 0 )
            {
                gmm->data[j] = (float2) (x[i], y[i]);
                ++j;
            }
        }

        // determine the number of clusters
        float bestValue = INFINITY;

        for ( char K = minClusters; K <= maxClusters; ++K )
        {
            // run each clustering sub-model
            bool success = GMM_fit(gmm, gmm->data, numSamples, K, gmm->labels);

            if ( !success )
            {
                continue;
            }

            // compute the criterion value of the sub-model
            float value = INFINITY;

            switch (criterion)
            {
            case AIC:
                value = GMM_computeAIC(K, 2, gmm->logL);
                break;
            case BIC:
                value = GMM_computeBIC(K, 2, gmm->logL, numSamples);
                break;
            case ICL:
                value = GMM_computeICL(K, 2, gmm->logL, numSamples, gmm->entropy);
                break;
            }

            // save the sub-model with the lowest criterion value
            if ( value < bestValue )
            {
                bestK = K;
                bestValue = value;

                // save labels for clean samples
                for ( int i = 0, j = 0; i < sampleSize; ++i )
                {
                    if ( labels[i] >= 0 )
                    {
                        labels[i] = gmm->labels[j];
                        ++j;
                    }
                }
            }
        }
    }

    // return number of clusters
    return bestK;
}
