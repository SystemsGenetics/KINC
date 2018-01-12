#include "genepair_gmm.h"



using namespace GenePair;





void GMM::Component::initialize(float pi, const Vector2& mu)
{
   // initialize pi and mu as given
   _pi = pi;
   _mu = mu;

   // Use identity covariance- assume dimensions are independent
   matrixInitIdentity(_sigma);

   // Initialize zero artifacts
   matrixInitZero(_sigmaInv);

   _normalizer = 0;
}






void GMM::Component::prepareCovariance()
{
   const int D = 2;

   // Compute inverse of Sigma once each iteration instead of
   // repeatedly for each calcLogMvNorm execution.
   float det;
   matrixInverse(_sigma, _sigmaInv, &det);

   if ( det <= 0 )
   {
      throw std::runtime_error("matrix inverse failed");
   }

   // Compute normalizer for multivariate normal distribution
   _normalizer = -0.5f * (D * log(2.0f * M_PI) + log(det));
}






void GMM::Component::calcLogMvNorm(const QVector<Vector2>& X, float *logP)
{
   // Here we are computing the probability density function of the multivariate
   // normal distribution conditioned on a single component for the set of points
   // given by X.
   //
   // P(x|k) = exp{ -0.5 * (x - mu)^T Sigma^{-} (x - mu) } / sqrt{ (2pi)^d det(Sigma) }

   for (int i = 0; i < X.size(); ++i)
   {
      // Let xm = (x - mu)
      Vector2 xm = X[i];
      vectorSubtract(xm, _mu);

      // Compute xm^T Sxm = xm^T S^-1 xm
      Vector2 Sxm;
      matrixProduct(_sigmaInv, xm, Sxm);

      float xmSxm = vectorDot(xm, Sxm);

      // Compute log(P) = normalizer - 0.5 * xm^T * S^-1 * xm
      logP[i] = _normalizer - 0.5f * xmSxm;
   }
}






void GMM::kmeans(const QVector<Vector2>& X)
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

      for (int i = 0; i < X.size(); ++i)
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






void GMM::calcLogMvNorm(const QVector<Vector2>& X, float *loggamma)
{
   const int N = X.size();
   const int K = _components.size();

   for ( int k = 0; k < K; ++k )
   {
      _components[k].calcLogMvNorm(X, &loggamma[k * N]);
   }
}






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






void GMM::performMStep(float *logpi, int K, float *loggamma, float *logGamma, float logGammaSum, const QVector<Vector2>& X)
{
   const int N = X.size();

   // update pi
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
      // Update mu
      Vector2& mu = _components[k]._mu;

      vectorInitZero(mu);

      for (int i = 0; i < N; ++i)
      {
         vectorAdd(mu, loggamma[k * N + i], X[i]);
      }

      vectorScale(mu, 1.0f / logGamma[k]);

      // Update sigma
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

      _components[k].prepareCovariance();
   }
}






QVector<cl_char> GMM::calcLabels(float *loggamma, int N, int K)
{
   QVector<cl_char> labels(N);

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

   return labels;
}






float GMM::calcEntropy(float *loggamma, int N)
{
   float E = 0;

   for ( int i = 0; i < N; ++i )
   {
      int k = _labels[i];

      E += log(loggamma[k * N + i]);
   }

   return E;
}






void GMM::fit(const QVector<Vector2>& X, int K)
{
   const int N = X.size();

   // initialize components
   _components.resize(K);

   for ( int k = 0; k < K; ++k )
   {
      // use uniform mixture proportion and randomly sampled mean
      int i = rand() % N;

      _components[k].initialize(1.0f / K, X[i]);
      _components[k].prepareCovariance();
   }

   // initialize means with k-means
   kmeans(X);

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

   try
   {
      for ( int t = 0; t < MAX_ITERATIONS; ++t )
      {
         // E step
         // compute gamma, log-likelihood
         calcLogMvNorm(X, loggamma);

         prevLogL = currentLogL;
         calcLogLikelihoodAndGammaNK(logpi, K, loggamma, N, &currentLogL);

         // check for convergence
         if ( fabs(currentLogL - prevLogL) < TOLERANCE )
         {
            break;
         }

         // M step
         // Let Gamma[k] = \Sum_i gamma[k, i]
         calcLogGammaK(loggamma, N, K, logGamma);

         float logGammaSum = calcLogGammaSum(logpi, K, logGamma);

         // Update parameters
         performMStep(logpi, K, loggamma, logGamma, logGammaSum, X);
      }

      // save outputs
      _success = true;
      _logL = currentLogL;
      _labels = calcLabels(loggamma, N, K);
      _entropy = calcEntropy(loggamma, N);
   }
   catch ( std::runtime_error& e )
   {
      _success = false;
   }

   delete[] logpi;
   delete[] loggamma;
   delete[] logGamma;
}
