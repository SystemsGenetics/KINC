#include "pairwise_gmm.h"



using namespace Pairwise;





/*!
 * Construct a Gaussian mixture model.
 *
 * @param emx
 */
GMM::GMM(ExpressionMatrix* emx):
   ClusteringModel(emx)
{
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

   // initialize precision to zero matrix
   matrixInitZero(_sigmaInv);

   // initialize normalizer term to 0
   _normalizer = 0;
}






/*!
 * Pre-compute the precision matrix and normalizer term for a mixture component.
 */
void GMM::Component::prepare()
{
   const int D = 2;

   // compute precision (inverse of covariance)
   float det;
   matrixInverse(_sigma, _sigmaInv, &det);

   if ( det <= 0 )
   {
      throw std::runtime_error("matrix inverse failed");
   }

   // compute normalizer term for multivariate normal distribution
   _normalizer = -0.5f * (D * log(2.0f * M_PI) + log(det));
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
   const float TOLERANCE = 1e-3;
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
         // determine the mean which is nearest to x_i
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
 * loggamma, the posterior probabilities for each component in the mixture model
 * and each sample in X, as well as logGamma, the sum of loggamma across samples,
 * and logGammaSum, another term used in the maximization step:
 *
 *   log(gamma_ki) = log(P(x_i|k)) - log(p(x))
 *
 *   log(p(x)) = a + log(sum(exp(logpi_k + logProb_ki) - a))
 *
 *   log(Gamma_k) = a + log(sum(exp(loggamma_ki) - a))
 *
 *   logGammaSum = a + log(sum(exp(logpi_k + logGamma_k) - a))
 *
 * @param X
 * @param N
 * @param logpi
 * @param loggamma
 * @param logGamma
 * @param logGammaSum
 */
float GMM::computeEStep(const QVector<Vector2>& X, int N, const float *logpi, float *loggamma, float *logGamma, float *logGammaSum)
{
   const int K = _components.size();

   // compute the log-probabilities for each component and each point in X
   float *logProb = loggamma;

   for ( int k = 0; k < K; ++k )
   {
      _components[k].computeLogProbNorm(X, N, &logProb[k * N]);
   }

   // compute loggamma and the log-likelihood
   float logL = 0.0;

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
      float sum = 0.0;
      for (int k = 0; k < K; ++k)
      {
         sum += exp(logpi[k] + logProb[k * N + i] - maxArg);
      }

      float logpx = maxArg + log(sum);

      // update loggamma_ki
      for (int k = 0; k < K; ++k)
      {
         loggamma[k * N + i] -= logpx;
      }

      // update log-likelihood
      logL += logpx;
   }

   // compute logGamma
   for (int k = 0; k < K; ++k)
   {
      const float *loggamma_k = &loggamma[k * N];

      // compute a = argmax(loggamma_ki, i)
      float maxArg = -INFINITY;
      for (int i = 0; i < N; ++i)
      {
         float arg = loggamma_k[i];
         if (maxArg < arg)
         {
            maxArg = arg;
         }
      }

      // compute logGamma_k
      float sum = 0;
      for (int i = 0; i < N; ++i)
      {
         sum += exp(loggamma_k[i] - maxArg);
      }

      logGamma[k] = maxArg + log(sum);
   }

   // compute logGammaSum
   float maxArg = -INFINITY;
   for (int k = 0; k < K; ++k)
   {
      float arg = logpi[k] + logGamma[k];
      if (maxArg < arg)
      {
         maxArg = arg;
      }
   }

   float sum = 0;
   for (int k = 0; k < K; ++k)
   {
      sum += exp(logpi[k] + logGamma[k] - maxArg);
   }

   *logGammaSum = maxArg + log(sum);

   // return log-likelihood
   return logL;
}






/*!
 * Perform the maximization step of the EM algorithm. In this step we update the
 * parameters of the the mixture model using loggamma, logGamma, and logGammaSum,
 * which are computed during the expectation step:
 *
 *   pi_k = exp(logpi_k + logGamma_k - logGammaSum)
 *
 *   mu_k = sum(exp(loggamma_ki) * x_i)) / exp(logGamma_k)
 *
 *   Sigma_k = sum(exp(loggamma_ki) * (x_i - mu_k) * (x_i - mu_k)^T) / exp(logGamma_k)
 *
 * @param X
 * @param N
 * @param logpi
 * @param loggamma
 * @param logGamma
 * @param logGammaSum
 */
void GMM::computeMStep(const QVector<Vector2>& X, int N, float *logpi, float *loggamma, float *logGamma, float logGammaSum)
{
   const int K = _components.size();

   // update mixture weight
   for (int k = 0; k < K; ++k)
   {
      logpi[k] += logGamma[k] - logGammaSum;

      _components[k]._pi = exp(logpi[k]);
   }

   // convert loggamma to gamma
   for (int k = 0; k < K; ++k)
   {
      for (int i = 0; i < N; ++i)
      {
         loggamma[k * N + i] = exp(loggamma[k * N + i]);
      }
   }

   // convert logGamma to Gamma
   for (int k = 0; k < K; ++k)
   {
      logGamma[k] = exp(logGamma[k]);
   }

   // update remaining parameters
   for (int k = 0; k < K; ++k)
   {
      // update mean
      Vector2& mu = _components[k]._mu;

      vectorInitZero(mu);

      for (int i = 0; i < N; ++i)
      {
         vectorAdd(mu, loggamma[k * N + i], X[i]);
      }

      vectorScale(mu, 1.0f / logGamma[k]);

      // update covariance matrix
      Matrix2x2& sigma = _components[k]._sigma;

      matrixInitZero(sigma);

      for (int i = 0; i < N; ++i)
      {
         // compute xm = (x_i - mu_k)
         Vector2 xm = X[i];
         vectorSubtract(xm, mu);

         // compute Sigma_ki = gamma_ki * (x_i - mu_k) (x_i - mu_k)^T
         Matrix2x2 outerProduct;
         matrixOuterProduct(xm, xm, outerProduct);

         matrixAdd(sigma, loggamma[k * N + i], outerProduct);
      }

      matrixScale(sigma, 1.0f / logGamma[k]);

      // pre-compute precision matrix and normalizer term
      _components[k].prepare();
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
 *   E = sum(sum(z_ki * log(gamma_ki))), z_ki = (y_i == k)
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

      E += log(gamma[k * N + i]);
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
   // initialize components
   _components.resize(K);

   for ( int k = 0; k < K; ++k )
   {
      // use uniform mixture weight and randomly sampled mean
      int i = rand() % N;

      _components[k].initialize(1.0f / K, X[i]);
      _components[k].prepare();
   }

   // initialize means with k-means
   initializeMeans(X, N);

   // initialize workspace
   float *logpi = new float[K];
   float *loggamma = new float[K * N];
   float *logGamma = new float[K];

   for (int k = 0; k < K; ++k)
   {
      logpi[k] = log(_components[k]._pi);
   }

   // run EM algorithm
   const int MAX_ITERATIONS = 100;
   const float TOLERANCE = 1e-8;
   float prevLogL = -INFINITY;
   float currLogL = -INFINITY;
   bool success;

   try
   {
      for ( int t = 0; t < MAX_ITERATIONS; ++t )
      {
         // perform E step
         float logGammaSum;

         prevLogL = currLogL;
         currLogL = computeEStep(X, N, logpi, loggamma, logGamma, &logGammaSum);

         // check for convergence
         if ( fabs(currLogL - prevLogL) < TOLERANCE )
         {
            break;
         }

         // perform M step
         computeMStep(X, N, logpi, loggamma, logGamma, logGammaSum);
      }

      // save outputs
      _logL = currLogL;
      computeLabels(loggamma, N, K, labels);
      _entropy = computeEntropy(loggamma, N, labels);

      success = true;
   }
   catch ( std::runtime_error& e )
   {
      success = false;
   }

   delete[] logpi;
   delete[] loggamma;
   delete[] logGamma;

   return success;
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

   return log(N) * p - 2 * logL;
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

   return log(N) * p - 2 * logL - 2 * E;
}
