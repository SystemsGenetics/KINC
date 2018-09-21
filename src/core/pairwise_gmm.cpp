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
 * P(x|k) = exp(-0.5 * (x - mu)^T Sigma^-1 (x - mu)) / sqrt((2pi)^d det(Sigma))
 *
 * Therefore the log-probability is:
 *
 * log(P(x|k)) = -0.5 * (x - mu)^T Sigma^-1 (x - mu) - 0.5 * (d * log(2pi) + log(det(Sigma)))
 *
 * @param X
 * @param N
 * @param logP
 */
void GMM::Component::calcLogMvNorm(const QVector<Vector2>& X, int N, float *logP)
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
void GMM::kmeans(const QVector<Vector2>& X, int N)
{
   const int K = _components.size();

   const int MAX_ITERATIONS = 20;
   const float TOLERANCE = 1e-3;
   float diff = 0;

   Vector2 MP[K];
   int counts[K];

   for (int t = 0; t < MAX_ITERATIONS && diff > TOLERANCE; ++t)
   {
      memset(MP, 0, K * sizeof(Vector2));
      memset(counts, 0, K * sizeof(int));

      for (int i = 0; i < N; ++i)
      {
         // arg min
         float minD = INFINITY;
         int minDk = 0;
         for (int k = 0; k < K; ++k)
         {
            float dist = vectorDiffNorm(X[i], _components[k]._mu);
            if (minD > dist)
            {
               minD = dist;
               minDk = k;
            }
         }

         vectorAdd(MP[minDk], X[i]);
         ++counts[minDk];
      }

      for (int k = 0; k < K; ++k)
      {
         vectorScale(MP[k], 1.0f / counts[k]);
      }

      diff = 0;
      for (int k = 0; k < K; ++k)
      {
         diff += vectorDiffNorm(MP[k], _components[k]._mu);
      }
      diff /= K;

      for (int k = 0; k < K; ++k)
      {
         _components[k]._mu = MP[k];
      }
   }
}






/*!
 * Compute the log of the pdf of a multivariate normal distribution for each
 * component in the mixture model and each point in X. The resulting matrix is
 * stored in loggamma.
 *
 * @param X
 * @param N
 * @param loggamma
 */
void GMM::calcLogMvNorm(const QVector<Vector2>& X, int N, float *loggamma)
{
   const int K = _components.size();

   for ( int k = 0; k < K; ++k )
   {
      _components[k].calcLogMvNorm(X, N, &loggamma[k * N]);
   }
}






/*!
 * Compute the log-likelihood of the mixture model as well as loggamma, the matrix
 * of log-probabilities for each component in the mixture model and each point in X.
 *
 * @param logpi
 * @param K
 * @param loggamma
 * @param N
 * @param logL
 */
void GMM::calcLogLikelihoodAndGammaNK(const float *logpi, int K, float *loggamma, int N, float *logL)
{
   *logL = 0.0;
   for (int i = 0; i < N; ++i)
   {
      float maxArg = -INFINITY;
      for (int k = 0; k < K; ++k)
      {
         const float logProbK = logpi[k] + loggamma[k * N + i];
         if (logProbK > maxArg)
         {
            maxArg = logProbK;
         }
      }

      float sum = 0.0;
      for (int k = 0; k < K; ++k)
      {
         const float logProbK = logpi[k] + loggamma[k * N + i];
         sum += exp(logProbK - maxArg);
      }

      const float logpx = maxArg + log(sum);
      *logL += logpx;
      for (int k = 0; k < K; ++k)
      {
         loggamma[k * N + i] += -logpx;
      }
   }
}






/*!
 * Compute logGamma.
 *
 * @param loggamma
 * @param N
 * @param K
 * @param logGamma
 */
void GMM::calcLogGammaK(const float *loggamma, int N, int K, float *logGamma)
{
   memset(logGamma, 0, K * sizeof(float));

   for (int k = 0; k < K; ++k)
   {
      const float *loggammak = &loggamma[k * N];

      float maxArg = -INFINITY;
      for (int i = 0; i < N; ++i)
      {
         const float loggammank = loggammak[i];
         if (loggammank > maxArg)
         {
            maxArg = loggammank;
         }
      }

      float sum = 0;
      for (int i = 0; i < N; ++i)
      {
         const float loggammank = loggammak[i];
         sum += exp(loggammank - maxArg);
      }

      logGamma[k] = maxArg + log(sum);
   }
}






/*!
 * Compute the sum of logGamma.
 *
 * @param logpi
 * @param K
 * @param logGamma
 */
float GMM::calcLogGammaSum(const float *logpi, int K, const float *logGamma)
{
   float maxArg = -INFINITY;
   for (int k = 0; k < K; ++k)
   {
      const float arg = logpi[k] + logGamma[k];
      if (arg > maxArg)
      {
         maxArg = arg;
      }
   }

   float sum = 0;
   for (int k = 0; k < K; ++k)
   {
      const float arg = logpi[k] + logGamma[k];
      sum += exp(arg - maxArg);
   }

   return maxArg + log(sum);
}






/*!
 * Update each component in the mixture model using loggamma and logGamma.
 *
 * @param logpi
 * @param K
 * @param loggamma
 * @param logGamma
 * @param logGammaSum
 * @param X
 * @param N
 */
void GMM::performMStep(float *logpi, int K, float *loggamma, float *logGamma, float logGammaSum, const QVector<Vector2>& X, int N)
{
   // update mixture weight
   for (int k = 0; k < K; ++k)
   {
      logpi[k] += logGamma[k] - logGammaSum;

      _components[k]._pi = exp(logpi[k]);
   }

   // convert loggamma / logGamma to gamma / Gamma to avoid duplicate exp(x) calls
   for (int k = 0; k < K; ++k)
   {
      for (int i = 0; i < N; ++i)
      {
         const int idx = k * N + i;
         loggamma[idx] = exp(loggamma[idx]);
      }
   }

   for (int k = 0; k < K; ++k)
   {
      logGamma[k] = exp(logGamma[k]);
   }

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
         // xm = (x - mu)
         Vector2 xm = X[i];
         vectorSubtract(xm, mu);

         // S_i = gamma_ik * (x - mu) (x - mu)^T
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
 * Compute the cluster labels of a dataset using loggamma.
 *
 * @param loggamma
 * @param N
 * @param K
 * @param labels
 */
void GMM::calcLabels(float *loggamma, int N, int K, QVector<qint8>& labels)
{
   for ( int i = 0; i < N; ++i )
   {
      int max_k = -1;
      float max_gamma = -INFINITY;

      for ( int k = 0; k < K; ++k )
      {
         if ( max_gamma < loggamma[k * N + i] )
         {
            max_k = k;
            max_gamma = loggamma[k * N + i];
         }
      }

      labels[i] = max_k;
   }
}






/*!
 * Compute the entropy between the mixture model and a dataset using loggamma
 * and the given cluster labels.
 *
 * @param loggamma
 * @param N
 * @param labels
 */
float GMM::calcEntropy(float *loggamma, int N, const QVector<qint8>& labels)
{
   float E = 0;

   for ( int i = 0; i < N; ++i )
   {
      int k = labels[i];

      E += log(loggamma[k * N + i]);
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
   kmeans(X, N);

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
   float currentLogL = -INFINITY;
   bool success;

   try
   {
      for ( int t = 0; t < MAX_ITERATIONS; ++t )
      {
         // E step
         // compute gamma, log-likelihood
         calcLogMvNorm(X, N, loggamma);

         prevLogL = currentLogL;
         calcLogLikelihoodAndGammaNK(logpi, K, loggamma, N, &currentLogL);

         // check for convergence
         if ( fabs(currentLogL - prevLogL) < TOLERANCE )
         {
            break;
         }

         // M step
         // compute Gamma[k] = \Sum_i gamma[k, i]
         calcLogGammaK(loggamma, N, K, logGamma);

         float logGammaSum = calcLogGammaSum(logpi, K, logGamma);

         // update parameters
         performMStep(logpi, K, loggamma, logGamma, logGammaSum, X, N);
      }

      // save outputs
      _logL = currentLogL;
      calcLabels(loggamma, N, K, labels);
      _entropy = calcEntropy(loggamma, N, labels);

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
