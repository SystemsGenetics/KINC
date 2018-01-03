#include <cassert>
#include <cfloat>
#include <gsl/gsl_cblas.h>

#include "genepair_gmm.h"



using namespace GenePair;





void choleskyDecomposition(const float *A, int D, float *L)
{
   // p. 157-158., Cholesky Factorization, 4.2 LU and Cholesky Factorizations,
   // Numerical Analysis by Kincaid, Cheney.

   // A is a real, symmetric, and positive definite D x D matrix
   assert(D > 0);
   assert(A != NULL);
   for (int i = 0; i < D; ++i) {
      for (int j = 0; j < D; ++j) {
         // Check that we are real valued
         float a = A[i * D + j];
         assert(a == a && fabsf(a) != INFINITY);

         // Check that we are symmetric
         float b = A[j * D + i];
         float absDiff = fabsf(a - b);
         assert(absDiff < 2.0 * FLT_EPSILON);
      }
   }

   // L is the resulting lower diagonal portion of A = LL^T
   assert(L != NULL);
   memset(L, 0, sizeof(float) * D * D);

   for (int k = 0; k < D; ++k) {
      float sum = 0;
      for (int s = 0; s < k; ++s) {
         const float l = L[k * D + s];
         const float ll = l * l;
         assert(ll == ll);
         assert(ll != -INFINITY);
         assert(ll != INFINITY);
         assert(ll >= 0);
         sum += ll;
      }

      assert(sum == sum);
      assert(sum != -INFINITY);
      assert(sum != INFINITY);
      assert(sum >= 0);

      sum = A[k * D + k] - sum;
      if (sum <= FLT_EPSILON) {
         // If this happens then we are not positive definite.
         throw std::runtime_error("A must be positive definite");
      }

      L[k * D + k] = sqrtf(sum);
      for (int i = k + 1; i < D; ++i) {
         float subsum = 0;
         for (int s = 0; s < k; ++s)
            subsum += L[i * D + s] * L[k * D + s];

         L[i * D + k] = (A[i * D + k] - subsum) / L[k * D + k];
      }
   }
}

void solvePositiveDefinite(const float *L, const float *B, float *X, int D, int N)
{
   // Want to solve the system given by: L(L^T)X = B where:
   //   L: D x D lower diagonal matrix
   //   X: D x N unknown
   //   B: D x N known
   //
   // Solve by first finding L Z = B, then L^T X = Z

   float *Z = (float *) calloc(N * D, sizeof(float));

   // 2015-09-23 GEL play the access of L into L(F)orward and L(B)ackward.
   // Found that sequential access improved runtime. 2017-03-24 GEL basically
   // pretend to carry out the forward and backward solvers, but to improve
   // runtime, load in L in sequential order ahead of time, so second time
   // around, we will have cached that data so the CPU will prefetch as needed.
   float *LF = (float *) malloc(D * D * sizeof(float));
   for (int i = 0, lf = 0; i < D; i++) {
      if (i > 0) {
         for (int j = 0; j <= i - 1; j++) {
            LF[lf++] = L[i * D + j];
         }
      }

      LF[lf++] = L[i * D + i];
   }

   float *LB = (float *) malloc(D * D * sizeof(float));
   for (int i = 0, lb = 0; i < D; ++i) {
      int ip = D - 1 - i;
      for (int j = ip + 1; j < D; j++) {
         LB[lb++] = L[j * D + ip];
      }

      LB[lb++] = L[ip * D + ip];
   }

   // Use forward subsitution to solve lower triangular matrix system Lz = b.
   // p. 150., Easy-to-Solve Systems, 4.2 LU and Cholesky Factorizations, Numerical Analysis by Kincaid, Cheney.
   for (int point = 0; point < N; ++point) {
      const float *b = &(B[point * D]);
      float *z = &(Z[point * D]);

      for (int i = 0, lf = 0; i < D; i++) {
         float sum = 0.0;
         if (i > 0) {
            for (int j = 0; j <= i - 1; j++) {
               sum += LF[lf++] * z[j];
            }
         }

         z[i] = (b[i] - sum) / LF[lf++];
      }
   }

   // use backward subsitution to solve L^T x = z
   // p. 150., Easy-to-Solve Systems, 4.2 LU and Cholesky Factorizations, Numerical Analysis by Kincaid, Cheney.
   for (int point = 0; point < N; ++point) {
      float *z = &(Z[point * D]);
      float *x = &(X[point * D]);

      for (int i = 0, lb = 0; i < D; i++) {
         int ip = D - 1 - i;

         float sum = 0;
         for (int j = ip + 1; j < D; j++)
            // Want A^T so switch switch i,j
            sum += LB[lb++] * x[j];

         x[ip] = (z[ip] - sum) / LB[lb++];
      }
   }

   free(LF);
   free(LB);
   free(Z);
}

float vecDiffNorm(const float *a, const float *b, int D)
{
   float dist = 0;
   for (int d = 0; d < D; ++d)
   {
      float distD = a[d] - b[d];
      distD *= distD;
      dist += distD;
   }
   return sqrtf(dist);
}






GMM::Component::Component()
{
   _D = 0;
   _pi = 0;
   _mu = nullptr;
   _sigma = nullptr;
   _sigmaL = nullptr;
   _normalizer = 0;
}






GMM::Component::~Component()
{
   free(_mu);
   free(_sigma);
   free(_sigmaL);
}






void GMM::Component::initialize(int D, float pi, float *mu)
{
   _D = D;
   _pi = pi;

   _mu = (float *) malloc(D * sizeof(float));
   memcpy(_mu, mu, D * sizeof(float));

   // Use identity covariance- assume dimensions are independent
   _sigma = (float *) calloc(D * D, sizeof(float));
   for ( int i = 0; i < D; ++i )
   {
      _sigma[i * D + i] = 1;
   }

   // Initialize zero artifacts
   _sigmaL = (float *) calloc(D * D, sizeof(float));
   _normalizer = 0;
}






void GMM::Component::prepareCovariance()
{
   // Perform cholesky factorization once each iteration instead of
   // repeadily for each normDist execution.
   choleskyDecomposition(_sigma, _D, _sigmaL);

   // log det(Sigma) = 2 sum log L_{i,i}
   float logDet = 1.0;
   for ( int i = 0; i < _D; ++i )
   {
      float diag = _sigmaL[i * _D + i];
      assert(diag > 0);
      logDet += logf(diag);
   }

   logDet *= 2;

   _normalizer = -0.5 * (_D * logf(2.0 * M_PI) + logDet);
}






void GMM::Component::logMvNormDist(const gsl_matrix_float *X, float *P)
{
   // Here we are computing the probability density function of the  multivariate
   // normal distribution conditioned on a single component for the set of points
   // given by X.
   //
   // P(x|component) = exp{ -0.5 * (x - mu)^T Sigma^{-} (x - mu) } / sqrt{ (2pi)^k det(Sigma) }
   //
   // Where Sigma and Mu are really Sigma_{component} Mu_{component}

   const int N = X->size1;
   const int D = X->size2;

   float *XM = (float *) malloc(N * D * sizeof(float));
   float *SXM = (float *) malloc(N * D * sizeof(float));
   float *innerProduct = (float *) malloc(N * sizeof(float));

   // Let XM = (x - m)
   for (int point = 0; point < N; ++point)
   {
      for (int dim = 0; dim < D; ++dim)
      {
         const int i = point * D + dim;
         XM[i] = X->data[i] - _mu[dim];
      }
   }

   // Sigma SXM = XM => Sigma^{-} XM = SXM
   solvePositiveDefinite(_sigmaL, XM, SXM, D, N);

   // XM^T SXM
   memset(innerProduct, 0, N * sizeof(float));
   for (int point = 0; point < N; ++point)
   {
      for (int dim = 0; dim < D; ++dim)
      {
         innerProduct[point] += XM[point * D + dim] * SXM[point * D + dim];
      }
   }

   // Compute P expf( -0.5 innerProduct ) / normalizer
   for (int point = 0; point < N; ++point)
   {
      // Normalizer already has negative sign on it.
      P[point] = -0.5 * innerProduct[point] + _normalizer;
      assert(P[point] == P[point]);
   }

   free(XM);
   free(SXM);
   free(innerProduct);
}






GMM::GMM()
{
   _K = 0;
   _success = false;
   _logL = INFINITY;
}





void GMM::initialize(const gsl_matrix_float *X, int K)
{
   // X is an N x D set of training data
   const int N = X->size1;
   const int D = X->size2;

   _K = K;
   _success = false;
   _logL = INFINITY;

   // initialize means randomly from X and with k-means
   float M[K * D];

   for ( int k = 0; k < K; ++k )
   {
      int i = rand() % N;

      for ( int j = 0; j < D; ++j )
      {
         M[k * D + j] = X->data[i * D + j];
      }
   }

   kmeans(X, M, K);

   // initialize components
   _components.resize(K);

   for ( int k = 0; k < K; ++k )
   {
      auto& component = _components[k];

      component.initialize(D, 1.0f / K, &M[k * D]);
      component.prepareCovariance();
   }
}






void GMM::kmeans(const gsl_matrix_float *X, float *M, int K)
{
   const int N = X->size1;
   const int D = X->size2;
   const float tolerance = 1e-3;
   float diff = 0;

   const int maxIterations = 20;

   float MP[K * D];
   int counts[K];

   for (int iteration = 0; iteration < maxIterations && diff > tolerance; ++iteration)
   {
      memset(MP, 0, K * D * sizeof(float));
      memset(counts, 0, K * sizeof(int));

      for (int i = 0; i < N; ++i)
      {
         const float *Xi = &X->data[i * D];

         // arg min
         float minD = INFINITY;
         int minDk = 0;
         for (int k = 0; k < K; ++k)
         {
            const float *Mk = &M[k * D];
            float dist = vecDiffNorm(Xi, Mk, D);
            if (minD > dist)
            {
               minD = dist;
               minDk = k;
            }
         }

         cblas_saxpy(D, 1.0f, Xi, 1, &M[minDk * D], 1);
         ++counts[minDk];
      }

      for (int k = 0; k < K; ++k)
      {
         cblas_sscal(D, 1.0f / counts[k], &MP[k * D], 1);
      }

      diff = 0;
      for (int k = 0; k < K; ++k)
      {
         diff += vecDiffNorm(&MP[k * D], &M[k * D], D);
      }
      diff /= (float) K;

      memcpy(M, MP, K * D * sizeof(float));
   }
}






void GMM::calcLogMvNorm(const gsl_matrix_float *X, float *logProb)
{
   const int N = X->size1;

   for ( int k = 0; k < _K; ++k )
   {
      auto& component = _components[k];

      component.logMvNormDist(X, &logProb[k * N]);
   }
}






void GMM::logLikelihoodAndGammaNK(const float *logpi, float *logProb, int N, float *logL)
{
   *logL = 0.0;
   for (int point = 0; point < N; ++point)
   {
      float maxArg = -INFINITY;
      for (int k = 0; k < _K; ++k)
      {
         const float logProbK = logpi[k] + logProb[k * N + point];
         if (logProbK > maxArg)
         {
            maxArg = logProbK;
         }
      }

      float sum = 0.0;
      for (int k = 0; k < _K; ++k)
      {
         const float logProbK = logpi[k] + logProb[k * N + point];
         sum += expf(logProbK - maxArg);
      }

      assert(sum >= 0);
      const float logpx = maxArg + logf(sum);
      *logL += logpx;
      for (int k = 0; k < _K; ++k)
      {
         logProb[k * N + point] += -logpx;
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
      for (int point = 0; point < N; ++point)
      {
         const float loggammank = loggammak[point];
         if (loggammank > maxArg)
         {
            maxArg = loggammank;
         }
      }

      float sum = 0;
      for (int point = 0; point < N; ++point)
      {
         const float loggammank = loggammak[point];
         sum += expf(loggammank - maxArg);
      }
      assert(sum >= 0);

      logGamma[k] = maxArg + logf(sum);
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
      sum += expf(arg - maxArg);
   }
   assert(sum >= 0);

   return maxArg + logf(sum);
}






void GMM::performMStep(float *logpi, float *loggamma, float *logGamma, const float logGammaSum, const gsl_matrix_float *X, float *outerProduct, float *xm)
{
   const int N = X->size1;
   const int D = X->size2;

   // update pi
   for (int k = 0; k < _K; ++k)
   {
      auto& component = _components[k];

      logpi[k] += logGamma[k] - logGammaSum;
      component._pi = expf(logpi[k]);
      assert(0 <= component._pi && component._pi <= 1);
   }

   // Convert loggamma and logGamma over to gamma and logGamma to avoid duplicate,
   //  and costly, expf(x) calls.
   for (int k = 0; k < _K; ++k)
   {
      for (int n = 0; n < N; ++n)
      {
         const int i = k * N + n;
         loggamma[i] = expf(loggamma[i]);
      }
   }

   for (int k = 0; k < _K; ++k)
   {
      logGamma[k] = expf(logGamma[k]);
   }

   // Update mu
   for (int k = 0; k < _K; ++k)
   {
      auto& component = _components[k];

      memset(component._mu, 0, D * sizeof(float));
      for (int point = 0; point < N; ++point)
      {
         for (int dim = 0; dim < D; ++dim)
         {
            component._mu[dim] += loggamma[k * N + point] * X->data[point * D + dim];
         }
      }

      for (int i = 0; i < D; ++i)
      {
         component._mu[i] /= logGamma[k];
      }
   }

   // Update sigma
   for (int k = 0; k < _K; ++k)
   {
      auto& component = _components[k];

      memset(component._sigma, 0, D * D * sizeof(float));
      for (int point = 0; point < N; ++point)
      {
         // (x - m)
         for (int dim = 0; dim < D; ++dim)
         {
            xm[dim] = X->data[point * D + dim] - component._mu[dim];
         }

         // (x - m) (x - m)^T
         for (int row = 0; row < D; ++row)
         {
            for (int column = 0; column < D; ++column)
            {
               outerProduct[row * D + column] = xm[row] * xm[column];
            }
         }

         for (int i = 0; i < D * D; ++i)
         {
            component._sigma[i] += loggamma[k * N + point] * outerProduct[i];
         }
      }

      for (int i = 0; i < D * D; ++i)
      {
         component._sigma[i] /= logGamma[k];
      }

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






void GMM::fit(const gsl_matrix_float *X, int K, int maxIterations)
{
   const int N = X->size1;
   const int D = X->size2;

   initialize(X, K);

   const float tolerance = 1e-8;
   int iteration = 0;
   float prevLogL = -INFINITY;
   float currentLogL = -INFINITY;

   float *logpi = (float *) calloc(K, sizeof(float));
   float *loggamma = (float *) calloc(N * K, sizeof(float));
   float *logGamma = (float *) calloc(K, sizeof(float));

   float *xm = (float *) calloc(D, sizeof(float));
   float *outerProduct = (float *) calloc(D * D, sizeof(float));

   for (int k = 0; k < K; ++k)
   {
      const float pik = _components[k]._pi;
      assert(pik >= 0);
      logpi[k] = logf(pik);
   }

   try
   {
      do
      {
         // --- E-Step ---

         // Compute gamma
         calcLogMvNorm(X, loggamma);

         prevLogL = currentLogL;
         logLikelihoodAndGammaNK(
            logpi,
            loggamma, N,
            &currentLogL
         );

         if ( fabsf(currentLogL - prevLogL) < tolerance )
         {
            break;
         }

         // Let Gamma[component] = \Sum_point gamma[component, point]
         calcLogGammaK(
            loggamma, N,
            logGamma
         );

         float logGammaSum = calcLogGammaSum(logpi, logGamma);

         // --- M-Step ---
         performMStep(
            logpi, loggamma, logGamma, logGammaSum,
            X,
            outerProduct, xm
         );
      } while ( ++iteration < maxIterations );

      // save outputs
      _success = true;
      _logL = currentLogL;
      _labels = calcLabels(loggamma, N);
   }
   catch ( std::runtime_error& e )
   {
      _success = false;
   }

   free(logpi);
   free(loggamma);
   free(logGamma);

   free(xm);
   free(outerProduct);
}
