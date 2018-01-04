#include <cassert>

#include "genepair_gmm.h"
#include "genepair_linalg.h"



using namespace GenePair;





GMM::Component::Component()
{
   _pi = 0;
   _mu = nullptr;
   _sigma = nullptr;
   _sigmaInv = nullptr;
   _normalizer = 0;
}






GMM::Component::~Component()
{
   delete[] _mu;
   delete[] _sigma;
   delete[] _sigmaInv;
}






void GMM::Component::initialize(float pi, float *mu)
{
   const int D = 2;

   // initialize pi and mu as given
   _pi = pi;
   _mu = mu;

   // Use identity covariance- assume dimensions are independent
   _sigma = new float[D * D];
   matrixInitIdentity(_sigma);

   // Initialize zero artifacts
   _sigmaInv = new float[D * D];
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

   // Compute normalizer for multivariate normal distribution
   _normalizer = -0.5f * (D * log(2.0f * M_PI) + log(det));
}






void GMM::Component::calcLogMvNorm(const float *X, int N, int D, float *logP)
{
   // Here we are computing the probability density function of the multivariate
   // normal distribution conditioned on a single component for the set of points
   // given by X.
   //
   // P(x|k) = exp{ -0.5 * (x - mu)^T Sigma^{-} (x - mu) } / sqrt{ (2pi)^d det(Sigma) }

   for (int i = 0; i < N; ++i)
   {
      // Let xm = (x - mu)
      float xm[D];
      vectorCopy(xm, &X[i * D]);
      vectorSubtract(xm, _mu);

      // Compute xm^T Sxm = xm^T S^-1 xm
      float Sxm[D];
      matrixProduct(_sigmaInv, xm, Sxm);

      float xmSxm = vectorDot(xm, Sxm);

      // Compute log(P) = normalizer - 0.5 * xm^T * S^-1 * xm
      logP[i] = _normalizer - 0.5f * xmSxm;
      assert(logP[i] == logP[i]);
   }
}






GMM::GMM()
{
   _K = 0;
   _success = false;
   _logL = INFINITY;
}





void GMM::initialize(const float *X, int N, int D, int K)
{
   // X is an N x D set of training data

   _K = K;
   _success = false;
   _logL = INFINITY;

   // initialize means randomly from X and with k-means
   QVector<float *> means(K);

   for ( int k = 0; k < K; ++k )
   {
      float *mu = new float[D];

      int i = rand() % N;
      vectorCopy(mu, &X[i * D]);

      means[k] = mu;
   }

   kmeans(X, N, D, means);

   // initialize components
   _components.resize(K);

   for ( int k = 0; k < K; ++k )
   {
      auto& component = _components[k];

      component.initialize(1.0f / K, means[k]);
      component.prepareCovariance();
   }
}






void GMM::kmeans(const float *X, int N, int D, const QVector<float *>& means)
{
   const int K = means.size();
   const float TOLERANCE = 1e-3;
   float diff = 0;

   const int MAX_ITERATIONS = 20;

   float MP[K * D];
   int counts[K];

   for (int t = 0; t < MAX_ITERATIONS && diff > TOLERANCE; ++t)
   {
      memset(MP, 0, K * D * sizeof(float));
      memset(counts, 0, K * sizeof(int));

      for (int i = 0; i < N; ++i)
      {
         const float *Xi = &X[i * D];

         // arg min
         float minD = INFINITY;
         int minDk = 0;
         for (int k = 0; k < K; ++k)
         {
            const float *Mk = means[k];
            float dist = vectorDiffNorm(Xi, Mk);
            if (minD > dist)
            {
               minD = dist;
               minDk = k;
            }
         }

         vectorAdd(means[minDk], Xi);
         ++counts[minDk];
      }

      for (int k = 0; k < K; ++k)
      {
         vectorScale(&MP[k * D], 1.0f / counts[k]);
      }

      diff = 0;
      for (int k = 0; k < K; ++k)
      {
         diff += vectorDiffNorm(&MP[k * D], means[k]);
      }
      diff /= K;

      for (int k = 0; k < K; ++k)
      {
         vectorCopy(means[k], &MP[k * D]);
      }
   }
}






void GMM::calcLogMvNorm(const float *X, int N, int D, float *logProb)
{
   for ( int k = 0; k < _K; ++k )
   {
      auto& component = _components[k];

      component.calcLogMvNorm(X, N, D, &logProb[k * N]);
   }
}






void GMM::calcLogLikelihoodAndGammaNK(const float *logpi, float *logProb, int N, float *logL)
{
   *logL = 0.0;
   for (int i = 0; i < N; ++i)
   {
      float maxArg = -INFINITY;
      for (int k = 0; k < _K; ++k)
      {
         const float logProbK = logpi[k] + logProb[k * N + i];
         if (logProbK > maxArg)
         {
            maxArg = logProbK;
         }
      }

      float sum = 0.0;
      for (int k = 0; k < _K; ++k)
      {
         const float logProbK = logpi[k] + logProb[k * N + i];
         sum += exp(logProbK - maxArg);
      }

      assert(sum >= 0);
      const float logpx = maxArg + log(sum);
      *logL += logpx;
      for (int k = 0; k < _K; ++k)
      {
         logProb[k * N + i] += -logpx;
      }
   }
}






void GMM::calcLogGammaK(const float *loggamma, int N, float *logGamma)
{
   memset(logGamma, 0, _K * sizeof(float));

   for (int k = 0; k < _K; ++k)
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
      assert(sum >= 0);

      logGamma[k] = maxArg + log(sum);
   }
}






float GMM::calcLogGammaSum(const float *logpi, const float *logGamma)
{
   float maxArg = -INFINITY;
   for (int k = 0; k < _K; ++k)
   {
      const float arg = logpi[k] + logGamma[k];
      if (arg > maxArg)
      {
         maxArg = arg;
      }
   }

   float sum = 0;
   for (int k = 0; k < _K; ++k)
   {
      const float arg = logpi[k] + logGamma[k];
      sum += exp(arg - maxArg);
   }
   assert(sum >= 0);

   return maxArg + log(sum);
}






void GMM::performMStep(float *logpi, float *loggamma, float *logGamma, float logGammaSum, const float *X, int N, int D)
{
   // update pi
   for (int k = 0; k < _K; ++k)
   {
      auto& component = _components[k];

      logpi[k] += logGamma[k] - logGammaSum;
      component._pi = exp(logpi[k]);
      assert(0 <= component._pi && component._pi <= 1);
   }

   // convert loggamma / logGamma to gamma / Gamma to avoid duplicate exp(x) calls
   for (int k = 0; k < _K; ++k)
   {
      for (int i = 0; i < N; ++i)
      {
         const int idx = k * N + i;
         loggamma[idx] = exp(loggamma[idx]);
      }
   }

   for (int k = 0; k < _K; ++k)
   {
      logGamma[k] = exp(logGamma[k]);
   }

   for (int k = 0; k < _K; ++k)
   {
      auto& component = _components[k];

      // Update mu
      float *mu = component._mu;

      vectorInitZero(mu);

      for (int i = 0; i < N; ++i)
      {
         vectorAdd(mu, loggamma[k * N + i], &X[i * D]);
      }

      vectorScale(mu, 1.0f / logGamma[k]);

      // Update sigma
      float *sigma = component._sigma;

      matrixInitZero(sigma);

      for (int i = 0; i < N; ++i)
      {
         // xm = (x - mu)
         float xm[D];
         vectorCopy(xm, &X[i * D]);
         vectorSubtract(xm, mu);

         // S_i = gamma_ik * (x - mu) (x - mu)^T
         float outerProduct[D * D];
         matrixOuterProduct(xm, xm, outerProduct);

         matrixAdd(sigma, loggamma[k * N + i], outerProduct);
      }

      matrixScale(sigma, 1.0f / logGamma[k]);

      component.prepareCovariance();
   }
}






QVector<int> GMM::calcLabels(float *loggamma, int N)
{
   QVector<int> labels(N);

   for ( int i = 0; i < N; i++ )
   {
      int max_j = -1;
      float max_gamma;

      for ( int j = 0; j < _K; j++ )
      {
         if ( max_j == -1 || max_gamma < loggamma[i * _K + j] )
         {
            max_j = j;
            max_gamma = loggamma[i * _K + j];
         }
      }

      labels[i] = max_j;
   }

   return labels;
}






void GMM::fit(const float *X, int N, int D, int K, int maxIterations)
{
   initialize(X, N, D, K);

   const float TOLERANCE = 1e-8;
   float prevLogL = -INFINITY;
   float currentLogL = -INFINITY;

   float *logpi = new float[K];
   float *loggamma = new float[N * K];
   float *logGamma = new float[K];

   for (int k = 0; k < K; ++k)
   {
      logpi[k] = log(_components[k]._pi);
   }

   try
   {
      for ( int t = 0; t < maxIterations; ++t )
      {
         // --- E-Step ---

         // Compute gamma
         calcLogMvNorm(X, N, D, loggamma);

         prevLogL = currentLogL;
         calcLogLikelihoodAndGammaNK(logpi, loggamma, N, &currentLogL);

         if ( fabs(currentLogL - prevLogL) < TOLERANCE )
         {
            break;
         }

         // Let Gamma[k] = \Sum_i gamma[k, i]
         calcLogGammaK(loggamma, N, logGamma);

         float logGammaSum = calcLogGammaSum(logpi, logGamma);

         // --- M-Step ---
         performMStep(logpi, loggamma, logGamma, logGammaSum, X, N, D);
      }

      // save outputs
      _success = true;
      _logL = currentLogL;
      _labels = calcLabels(loggamma, N);
   }
   catch ( std::runtime_error& e )
   {
      _success = false;
   }

   delete[] logpi;
   delete[] loggamma;
   delete[] logGamma;
}
