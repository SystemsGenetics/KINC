#include "genepair_kmeans.h"
#include "genepair_linalg.h"



using namespace GenePair;






KMeans::~KMeans()
{
   for ( gsl_vector_float *mu : _means ) {
      gsl_vector_float_free(mu);
   }
}






void KMeans::initialize(const gsl_matrix_float *X, int K)
{
   // initialize number of clusters
   _K = K;

   // initialize means randomly from X
   const int N = X->size1;
   const int D = X->size2;

   for ( int k = 0; k < K; ++k )
   {
      gsl_vector_float *mu = gsl_vector_float_alloc(D);

      int i = rand() % N;
      vectorCopy(mu->data, &X->data[i * D]);

      _means.push_back(mu);
   }
}






void KMeans::fit(const gsl_matrix_float *X, int K)
{
   // initialize model
   initialize(X, K);

   // iterate K means until convergence
   const int N = X->size1;
   const int D = X->size2;

   QVector<int> y;
   QVector<int> y_next(N);

   while ( true )
   {
      // E step
      for ( int i = 0; i < N; ++i )
      {
         // find k that minimizes norm(x_i - mu_k)
         int min_k = -1;
         float min_dist;

         for ( int k = 0; k < K; ++k )
         {
            const float *x_i = &X->data[i * D];
            const float *mu_k = _means[k]->data;

            float dist = vectorDiffNorm(x_i, mu_k);

            if ( min_k == -1 || dist < min_dist )
            {
               min_k = k;
               min_dist = dist;
            }
         }

         y_next[i] = min_k;
      }

      // check for convergence
      if ( y == y_next )
      {
         break;
      }

      y = y_next;

      // M step
      for ( int k = 0; k < K; ++k )
      {
         // compute mu_k = mean of all x_i in cluster k
         float *mu_k = _means[k]->data;
         int n_k = 0;

         vectorInitZero(mu_k);

         for ( int i = 0; i < N; ++i )
         {
            const float *x_i = &X->data[i * D];

            if ( y[i] == k )
            {
               vectorAdd(mu_k, x_i);
               n_k++;
            }
         }

         vectorScale(mu_k, 1.0f / n_k);
      }
   }

   // save outputs
   _logL = computeLogLikelihood(X);
   _labels = y;
}






float KMeans::computeLogLikelihood(const gsl_matrix_float *X)
{
   const int N = X->size1;
   const int D = X->size2;

   // compute within-class scatter
   float S = 0;

   for ( int k = 0; k < _K; ++k )
   {
      for ( int i = 0; i < N; ++i )
      {
         if ( _labels[i] != k )
         {
            continue;
         }

         const float *mu_k = _means[k]->data;
         const float *x_i = &X->data[i * D];

         float dist = vectorDiffNorm(x_i, mu_k);

         S += dist * dist;
      }
   }

   return -S;
}
