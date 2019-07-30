
// #include "fetchpair.cl"
// #include "linalg.cl"
// #include "outlier.cl"






typedef struct
{
   float pi;
   Vector2 mu;
   Matrix2x2 sigma;
   Matrix2x2 sigmaInv;
   float normalizer;
} Component;






typedef struct
{
   __global Component *components;
   int K;
   float logL;
   float entropy;
   __global Vector2 *_Mu;
   __global int *_counts;
   __global float *_logpi;
   __global float *_gamma;
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
 * Initialize a mixture component with the given mixture weight and mean.
 *
 * @param component
 * @param pi
 * @param mu
 */
void GMM_Component_initialize(
   __global Component *component,
   float pi,
   __global const Vector2 *mu)
{
   // initialize mixture weight and mean
   component->pi = pi;
   component->mu = *mu;

   // initialize covariance to identity matrix
   matrixInitIdentity(&component->sigma);
}






/*!
 * Pre-compute the precision matrix and normalizer term for a mixture component.
 *
 * @param component
 */
bool GMM_Component_prepare(__global Component *component)
{
   const int D = 2;

   // compute precision (inverse of covariance)
   float det;
   matrixInverse(&component->sigma, &component->sigmaInv, &det);

   // compute normalizer term for multivariate normal distribution
   component->normalizer = -0.5f * (D * log(2.0f * M_PI) + log(det));

   // return failure if matrix inverse failed
   return !(det <= 0 || isnan(det));
}






/*!
 * Compute the log of the probability density function of the multivariate normal
 * distribution conditioned on a single component for each point in X:
 *
 *   P(x|k) = exp(-0.5 * (x - mu)^T Sigma^-1 (x - mu)) / sqrt((2pi)^d det(Sigma))
 *
 * Therefore the log-probability is:
 *
 *   log(P(x|k)) = -0.5 * (x - mu)^T Sigma^-1 (x - mu) - 0.5 * (d * log(2pi) + log(det(Sigma)))
 *
 * @param component
 * @param X
 * @param N
 * @param logP
 */
void GMM_Component_computeLogProbNorm(
   __global const Component *component,
   __global const Vector2 *X, int N,
   __global float *logP)
{
   for ( int i = 0; i < N; ++i )
   {
      // compute xm = (x - mu)
      Vector2 xm = X[i];
      vectorSubtract(&xm, &component->mu);

      // compute Sxm = Sigma^-1 xm
      Vector2 Sxm;
      matrixProduct(&component->sigmaInv, &xm, &Sxm);

      // compute xmSxm = xm^T Sigma^-1 xm
      float xmSxm = vectorDot(&xm, &Sxm);

      // compute log(P) = normalizer - 0.5 * xm^T * Sigma^-1 * xm
      logP[i] = component->normalizer - 0.5f * xmSxm;
   }
}






/*!
 * Initialize the mean of each component in the mixture model using k-means
 * clustering.
 *
 * @param gmm
 * @param X
 * @param N
 */
void GMM_initializeMeans(GMM *gmm, __global const Vector2 *X, int N)
{
   const int K = gmm->K;

   const int MAX_ITERATIONS = 20;
   const float TOLERANCE = 1e-3f;
   float diff = 0;

   // initialize workspace
   __global Vector2 *Mu = gmm->_Mu;
   __global int *counts = gmm->_counts;

   for ( int t = 0; t < MAX_ITERATIONS && diff > TOLERANCE; ++t )
   {
      // compute mean and sample count for each component
      for ( int k = 0; k < K; ++k )
      {
         vectorInitZero(&Mu[k]);
         counts[k] = 0;
      }

      for ( int i = 0; i < N; ++i )
      {
         // determine the component mean which is nearest to x_i
         float min_dist = INFINITY;
         int min_k = 0;
         for ( int k = 0; k < K; ++k )
         {
            float dist = vectorDiffNorm(&X[i], &gmm->components[k].mu);
            if ( min_dist > dist )
            {
               min_dist = dist;
               min_k = k;
            }
         }

         // update mean and sample count
         vectorAdd(&Mu[min_k], &X[i]);
         ++counts[min_k];
      }

      // scale each mean by its sample count
      for ( int k = 0; k < K; ++k )
      {
         vectorScale(&Mu[k], 1.0f / counts[k]);
      }

      // compute the total change of all means
      diff = 0;
      for ( int k = 0; k < K; ++k )
      {
         diff += vectorDiffNorm(&Mu[k], &gmm->components[k].mu);
      }
      diff /= K;

      // update component means
      for ( int k = 0; k < K; ++k )
      {
         gmm->components[k].mu = Mu[k];
      }
   }
}






/*!
 * Perform the expectation step of the EM algorithm. In this step we compute
 * gamma, the posterior probabilities for each component in the mixture model
 * and each sample in X, as well as the log-likelihood of the model:
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
 */
float GMM_computeEStep(GMM *gmm, __global const Vector2 *X, int N)
{
   const int K = gmm->K;

   // compute logpi
   for ( int k = 0; k < K; ++k )
   {
      gmm->_logpi[k] = log(gmm->components[k].pi);
   }

   // compute the log-probability for each component and each point in X
   __global float *logProb = gmm->_gamma;

   for ( int k = 0; k < K; ++k )
   {
      GMM_Component_computeLogProbNorm(&gmm->components[k], X, N, &logProb[k * N]);
   }

   // compute gamma and log-likelihood
   float logL = 0;

   for ( int i = 0; i < N; ++i )
   {
      // compute a = argmax(logpi_k + logProb_ki, k)
      float maxArg = -INFINITY;
      for ( int k = 0; k < K; ++k )
      {
         float arg = gmm->_logpi[k] + logProb[k * N + i];
         if ( maxArg < arg )
         {
            maxArg = arg;
         }
      }

      // compute logpx
      float sum = 0;
      for ( int k = 0; k < K; ++k )
      {
         sum += exp(gmm->_logpi[k] + logProb[k * N + i] - maxArg);
      }

      float logpx = maxArg + log(sum);

      // compute gamma_ki
      for ( int k = 0; k < K; ++k )
      {
         gmm->_gamma[k * N + i] += gmm->_logpi[k] - logpx;
         gmm->_gamma[k * N + i] = exp(gmm->_gamma[k * N + i]);
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
 */
void GMM_computeMStep(GMM *gmm, __global const Vector2 *X, int N)
{
   const int K = gmm->K;

   for ( int k = 0; k < K; ++k )
   {
      // compute n_k = sum(gamma_ki)
      float n_k = 0;

      for ( int i = 0; i < N; ++i )
      {
         n_k += gmm->_gamma[k * N + i];
      }

      // update mixture weight
      gmm->components[k].pi = n_k / N;

      // update mean
      __global Vector2 *mu = &gmm->components[k].mu;

      vectorInitZero(mu);

      for ( int i = 0; i < N; ++i )
      {
         vectorAddScaled(mu, gmm->_gamma[k * N + i], &X[i]);
      }

      vectorScale(mu, 1.0f / n_k);

      // update covariance matrix
      __global Matrix2x2 *sigma = &gmm->components[k].sigma;

      matrixInitZero(sigma);

      for ( int i = 0; i < N; ++i )
      {
         // compute xm = (x_i - mu_k)
         Vector2 xm = X[i];
         vectorSubtract(&xm, mu);

         // compute Sigma_ki = gamma_ki * (x_i - mu_k) (x_i - mu_k)^T
         matrixAddOuterProduct(sigma, gmm->_gamma[k * N + i], &xm);
      }

      matrixScale(sigma, 1.0f / n_k);
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
   __global const float *gamma, int N, int K,
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
   __global const float *gamma, int N,
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
   __global const Vector2 *X, int N, int K,
   __global char *labels)
{
   // initialize random state
   ulong state = 1;

   // initialize components
   gmm->K = K;

   for ( int k = 0; k < K; ++k )
   {
      // use uniform mixture weight and randomly sampled mean
      int i = myrand(&state) % N;

      GMM_Component_initialize(&gmm->components[k], 1.0f / K, &X[i]);
   }

   // initialize means with k-means
   GMM_initializeMeans(gmm, X, N);

   // run EM algorithm
   const int MAX_ITERATIONS = 100;
   const float TOLERANCE = 1e-8f;
   float prevLogL = -INFINITY;
   float currLogL = -INFINITY;

   for ( int t = 0; t < MAX_ITERATIONS; ++t )
   {
      // pre-compute precision matrix and normalizer term
      bool success = true;

      for ( int k = 0; k < K; k++ )
      {
         success = success && GMM_Component_prepare(&gmm->components[k]);
      }

      // return failure if matrix inverse failed
      if ( !success )
      {
         return false;
      }

      // perform E step
      prevLogL = currLogL;
      currLogL = GMM_computeEStep(gmm, X, N);

      // check for convergence
      if ( fabs(currLogL - prevLogL) < TOLERANCE )
      {
         break;
      }

      // perform M step
      GMM_computeMStep(gmm, X, N);
   }

   // save outputs
   gmm->logL = currLogL;
   GMM_computeLabels(gmm->_gamma, N, K, labels);
   gmm->entropy = GMM_computeEntropy(gmm->_gamma, N, labels);

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
 * @param numPairs
 * @param expressions
 * @param sampleSize
 * @param in_index
 * @param minSamples
 * @param minClusters
 * @param maxClusters
 * @param criterion
 * @param out_K
 * @param out_labels
 */
__kernel void GMM_compute(
   int numPairs,
   __global const float *expressions,
   int sampleSize,
   __global const int2 *in_index,
   int minSamples,
   char minClusters,
   char maxClusters,
   Criterion criterion,
   __global Vector2 *work_X,
   __global int *work_N,
   __global char *work_labels,
   __global Component *work_components,
   __global Vector2 *work_MP,
   __global int *work_counts,
   __global float *work_logpi,
   __global float *work_gamma,
   __global char *out_K,
   __global char *out_labels)
{
   int i = get_global_id(0);

   if ( i >= numPairs )
   {
      return;
   }

   // initialize workspace variables
   int2 index = in_index[i];
   __global const float *x = &expressions[index.x * sampleSize];
   __global const float *y = &expressions[index.y * sampleSize];
   __global Vector2 *X = &work_X[i * sampleSize];
   int numSamples = work_N[i];
   __global char *labels = &work_labels[i * sampleSize];
   __global Component *components = &work_components[i * maxClusters];
   __global Vector2 *Mu = &work_MP[i * maxClusters];
   __global int *counts = &work_counts[i * maxClusters];
   __global float *logpi = &work_logpi[i * maxClusters];
   __global float *gamma = &work_gamma[i * maxClusters * sampleSize];
   __global char *bestK = &out_K[i];
   __global char *bestLabels = &out_labels[i * sampleSize];

   // initialize GMM struct
   GMM gmm = {
      .components = components,
      ._Mu = Mu,
      ._counts = counts,
      ._logpi = logpi,
      ._gamma = gamma
   };

   // perform clustering only if there are enough samples
   *bestK = 0;

   if ( numSamples >= minSamples )
   {
      // extract clean samples from data array
      for ( int i = 0, j = 0; i < sampleSize; ++i )
      {
         if ( bestLabels[i] >= 0 )
         {
            X[j] = (float2) (x[i], y[i]);
            ++j;
         }
      }

      // determine the number of clusters
      float bestValue = INFINITY;

      for ( char K = minClusters; K <= maxClusters; ++K )
      {
         // run each clustering sub-model
         bool success = GMM_fit(&gmm, X, numSamples, K, labels);

         if ( !success )
         {
            continue;
         }

         // compute the criterion value of the sub-model
         float value = INFINITY;

         switch (criterion)
         {
         case AIC:
            value = GMM_computeAIC(K, 2, gmm.logL);
            break;
         case BIC:
            value = GMM_computeBIC(K, 2, gmm.logL, numSamples);
            break;
         case ICL:
            value = GMM_computeICL(K, 2, gmm.logL, numSamples, gmm.entropy);
            break;
         }

         // save the sub-model with the lowest criterion value
         if ( value < bestValue )
         {
            *bestK = K;
            bestValue = value;

            // save labels for clean samples
            for ( int i = 0, j = 0; i < sampleSize; ++i )
            {
               if ( bestLabels[i] >= 0 )
               {
                  bestLabels[i] = labels[j];
                  ++j;
               }
            }
         }
      }
   }
}
