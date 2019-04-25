#include "pairwise_gmm.h"



using namespace Pairwise;





/*!
 * Implementation of rand(), taken from POSIX example.
 *
 * @param state
 */
int myrand(unsigned long *state)
{
   *state = (*state) * 1103515245 + 12345;
   return ((unsigned)((*state)/65536) % 32768);
}






/*!
 * Construct a Gaussian mixture model.
 *
 * @param emx
 * @param maxClusters
 */
GMM::GMM(ExpressionMatrix* emx, qint8 maxClusters)
{
   // pre-allocate workspace
   _data.resize(emx->sampleSize());
   _labels.resize(emx->sampleSize());
   _components.reserve(maxClusters);
   _gamma = new float[maxClusters * emx->sampleSize()];
}






/*!
 * Destruct a Gaussian mixture model.
 */
GMM::~GMM()
{
   delete[] _gamma;
}






/*!
 * Initialize a mixture component with the given mixture weight and mean.
 *
 * @param pi
 * @param mu
 */
void GMM::Component::initialize(float pi, const Vector2& mu)
{
   // initialize mixture weight and mean
   _pi = pi;
   _mu = mu;

   // initialize covariance to identity matrix
   matrixInitIdentity(_sigma);
}






/*!
 * Pre-compute the precision matrix and normalizer term for a mixture component.
 */
bool GMM::Component::prepare()
{
   const int D = 2;

   // compute precision (inverse of covariance)
   float det;
   matrixInverse(_sigma, _sigmaInv, &det);

   // compute normalizer term for multivariate normal distribution
   _normalizer = -0.5f * (D * logf(2.0f * M_PI) + logf(det));

   // return failure if matrix inverse failed
   return !(det <= 0);
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
 * @param X
 * @param N
 * @param logP
 */
void GMM::Component::computeLogProbNorm(const QVector<Vector2>& X, int N, float *logP)
{
   for (int i = 0; i < N; ++i)
   {
      // compute xm = (x - mu)
      Vector2 xm = X[i];
      vectorSubtract(xm, _mu);

      // compute Sxm = Sigma^-1 xm
      Vector2 Sxm;
      matrixProduct(_sigmaInv, xm, Sxm);

      // compute xmSxm = xm^T Sigma^-1 xm
      float xmSxm = vectorDot(xm, Sxm);

      // compute log(P) = normalizer - 0.5 * xm^T * Sigma^-1 * xm
      logP[i] = _normalizer - 0.5f * xmSxm;
   }
}






/*!
 * Initialize the mean of each component in the mixture model using k-means
 * clustering.
 *
 * @param X
 * @param N
 */
void GMM::initializeMeans(const QVector<Vector2>& X, int N)
{
   const int K = _components.size();

   const int MAX_ITERATIONS = 20;
   const float TOLERANCE = 1e-3f;
   float diff = 0;

   // initialize workspace
   Vector2 Mu[K];
   int counts[K];

   for (int t = 0; t < MAX_ITERATIONS && diff > TOLERANCE; ++t)
   {
      // compute mean and sample count for each component
      memset(Mu, 0, K * sizeof(Vector2));
      memset(counts, 0, K * sizeof(int));

      for (int i = 0; i < N; ++i)
      {
         // determine the component mean which is nearest to x_i
         float min_dist = INFINITY;
         int min_k = 0;
         for (int k = 0; k < K; ++k)
         {
            float dist = vectorDiffNorm(X[i], _components[k]._mu);
            if (min_dist > dist)
            {
               min_dist = dist;
               min_k = k;
            }
         }

         // update mean and sample count
         vectorAdd(Mu[min_k], X[i]);
         ++counts[min_k];
      }

      // scale each mean by its sample count
      for (int k = 0; k < K; ++k)
      {
         vectorScale(Mu[k], 1.0f / counts[k]);
      }

      // compute the total change of all means
      diff = 0;
      for (int k = 0; k < K; ++k)
      {
         diff += vectorDiffNorm(Mu[k], _components[k]._mu);
      }
      diff /= K;

      // update component means
      for (int k = 0; k < K; ++k)
      {
         _components[k]._mu = Mu[k];
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
 * @param X
 * @param N
 */
float GMM::computeEStep(const QVector<Vector2>& X, int N)
{
   const int K = _components.size();

   // compute logpi
   float logpi[K];

   for (int k = 0; k < K; ++k)
   {
      logpi[k] = logf(_components[k]._pi);
   }

   // compute the log-probability for each component and each point in X
   float *logProb = _gamma;

   for ( int k = 0; k < K; ++k )
   {
      _components[k].computeLogProbNorm(X, N, &logProb[k * N]);
   }

   // compute gamma and log-likelihood
   float logL = 0;

   for (int i = 0; i < N; ++i)
   {
      // compute a = argmax(logpi_k + logProb_ki, k)
      float maxArg = -INFINITY;
      for (int k = 0; k < K; ++k)
      {
         float arg = logpi[k] + logProb[k * N + i];
         if (maxArg < arg)
         {
            maxArg = arg;
         }
      }

      // compute logpx
      float sum = 0;
      for (int k = 0; k < K; ++k)
      {
         sum += expf(logpi[k] + logProb[k * N + i] - maxArg);
      }

      float logpx = maxArg + logf(sum);

      // compute gamma_ki
      for (int k = 0; k < K; ++k)
      {
         _gamma[k * N + i] += logpi[k] - logpx;
         _gamma[k * N + i] = expf(_gamma[k * N + i]);
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
 * @param X
 * @param N
 */
void GMM::computeMStep(const QVector<Vector2>& X, int N)
{
   const int K = _components.size();

   for (int k = 0; k < K; ++k)
   {
      // compute n_k = sum(gamma_ki)
      float n_k = 0;

      for (int i = 0; i < N; ++i)
      {
         n_k += _gamma[k * N + i];
      }

      // update mixture weight
      _components[k]._pi = n_k / N;

      // update mean
      Vector2& mu = _components[k]._mu;

      vectorInitZero(mu);

      for (int i = 0; i < N; ++i)
      {
         vectorAdd(mu, _gamma[k * N + i], X[i]);
      }

      vectorScale(mu, 1.0f / n_k);

      // update covariance matrix
      Matrix2x2& sigma = _components[k]._sigma;

      matrixInitZero(sigma);

      for (int i = 0; i < N; ++i)
      {
         // compute xm = (x_i - mu_k)
         Vector2 xm = X[i];
         vectorSubtract(xm, mu);

         // compute Sigma_ki = gamma_ki * (x_i - mu_k) (x_i - mu_k)^T
         matrixAddOuterProduct(sigma, _gamma[k * N + i], xm);
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
void GMM::computeLabels(const float *gamma, int N, int K, QVector<qint8>& labels)
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
float GMM::computeEntropy(const float *gamma, int N, const QVector<qint8>& labels)
{
   float E = 0;

   for ( int i = 0; i < N; ++i )
   {
      int k = labels[i];

      E -= logf(gamma[k * N + i]);
   }

   return E;
}






/*!
 * Fit the mixture model to a pairwise data array and compute the output cluster
 * labels for the data. The data array should only contain clean samples.
 *
 * @param X
 * @param N
 * @param K
 * @param labels
 */
bool GMM::fit(const QVector<Vector2>& X, int N, int K, QVector<qint8>& labels)
{
   // initialize random state
   unsigned long state = 1;

   // initialize components
   _components.resize(K);

   for ( int k = 0; k < K; ++k )
   {
      // use uniform mixture weight and randomly sampled mean
      int i = myrand(&state) % N;

      _components[k].initialize(1.0f / K, X[i]);
   }

   // initialize means with k-means
   initializeMeans(X, N);

   // run EM algorithm
   const int MAX_ITERATIONS = 100;
   const float TOLERANCE = 1e-8f;
   float prevLogL = -INFINITY;
   float currLogL = -INFINITY;

   for ( int t = 0; t < MAX_ITERATIONS; ++t )
   {
      // pre-compute precision matrix and normalizer term
      bool success = true;

      for ( int k = 0; k < K; ++k )
      {
         success = success && _components[k].prepare();
      }

      // return failure if matrix inverse failed
      if ( !success )
      {
         return false;
      }

      // perform E step
      prevLogL = currLogL;
      currLogL = computeEStep(X, N);

      // check for convergence
      if ( fabs(currLogL - prevLogL) < TOLERANCE )
      {
         break;
      }

      // perform M step
      computeMStep(X, N);
   }

   // save outputs
   _logL = currLogL;
   computeLabels(_gamma, N, K, labels);
   _entropy = computeEntropy(_gamma, N, labels);

   return true;
}






/*!
 * Compute the Akaike Information Criterion of a Gaussian mixture model.
 *
 * @param K
 * @param D
 * @param logL
 */
float GMM::computeAIC(int K, int D, float logL)
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
float GMM::computeBIC(int K, int D, float logL, int N)
{
   int p = K * (1 + D + D * D);

   return logf(N) * p - 2 * logL;
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
float GMM::computeICL(int K, int D, float logL, int N, float E)
{
   int p = K * (1 + D + D * D);

   return logf(N) * p - 2 * logL + 2 * E;
}






/*!
 * Determine the number of clusters in a pairwise data array. Several sub-models,
 * each one having a different number of clusters, are fit to the data and the
 * sub-model with the best criterion value is selected.
 *
 * @param expressions
 * @param index
 * @param numSamples
 * @param labels
 * @param minSamples
 * @param minClusters
 * @param maxClusters
 * @param criterion
 */
qint8 GMM::compute(
   const std::vector<float>& expressions,
   const Index& index,
   int numSamples,
   QVector<qint8>& labels,
   int minSamples,
   qint8 minClusters,
   qint8 maxClusters,
   Criterion criterion)
{
   // index into gene expressions
   const float *x = &expressions[index.getX() * labels.size()];
   const float *y = &expressions[index.getY() * labels.size()];

   // perform clustering only if there are enough samples
   qint8 bestK = 0;

   if ( numSamples >= minSamples )
   {
      // extract clean samples from data array
      for ( int i = 0, j = 0; i < labels.size(); ++i )
      {
         if ( labels[i] >= 0 )
         {
            _data[j] = { x[i], y[i] };
            ++j;
         }
      }

      // determine the number of clusters
      float bestValue = INFINITY;

      for ( qint8 K = minClusters; K <= maxClusters; ++K )
      {
         // run each clustering sub-model
         bool success = fit(_data, numSamples, K, _labels);

         if ( !success )
         {
            continue;
         }

         // compute the criterion value of the sub-model
         float value = INFINITY;

         switch (criterion)
         {
         case Criterion::AIC:
            value = computeAIC(K, 2, _logL);
            break;
         case Criterion::BIC:
            value = computeBIC(K, 2, _logL, numSamples);
            break;
         case Criterion::ICL:
            value = computeICL(K, 2, _logL, numSamples, _entropy);
            break;
         }

         // save the sub-model with the lowest criterion value
         if ( value < bestValue )
         {
            bestK = K;
            bestValue = value;

            // save labels for clean samples
            for ( int i = 0, j = 0; i < labels.size(); ++i )
            {
               if ( labels[i] >= 0 )
               {
                  labels[i] = _labels[j];
                  ++j;
               }
            }
         }
      }
   }

   return bestK;
}
