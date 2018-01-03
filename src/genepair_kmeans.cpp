#include <gsl/gsl_cblas.h>

#include "genepair_kmeans.h"



using namespace GenePair;






float vector_diff_norm(const float *a, const float *b, int n)
{
   float dist = 0;
   for ( int i = 0; i < n; ++i )
   {
      float diff = a[i] - b[i];
      dist += diff * diff;
   }

   return sqrt(dist);
}






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
      memcpy(mu->data, &X->data[i * D], D * sizeof(float));

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

            float dist = vector_diff_norm(x_i, mu_k, D);

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

         memset(mu_k, 0, D * sizeof(float));

         for ( int i = 0; i < N; ++i )
         {
            const float *x_i = &X->data[i * D];

            if ( y[i] == k )
            {
               cblas_saxpy(D, 1.0f, x_i, 1, mu_k, 1);
               n_k++;
            }
         }

         cblas_sscal(D, 1.0f / n_k, mu_k, 1);
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

         float dist = vector_diff_norm(x_i, mu_k, D);

         S += dist * dist;
      }
   }

   return -S;
}
