#include "genepair_kmeans.h"



using namespace GenePair;






bool KMeans::fit(const QVector<Vector2>& X, int N, int K, QVector<qint8>& labels)
{
   const int NUM_INITS = 10;
   const int MAX_ITERATIONS = 300;

   // repeat with several initializations
   _logL = -INFINITY;

   for ( int init = 0; init < NUM_INITS; ++init )
   {
      // initialize means randomly from X
      _means.resize(K);

      for ( int k = 0; k < K; ++k )
      {
         int i = rand() % N;
         _means[k] = X[i];
      }

      // iterate K means until convergence
      QVector<qint8> y(N);
      QVector<qint8> y_next(N);

      for ( int t = 0; t < MAX_ITERATIONS; ++t )
      {
         // compute new labels
         for ( int i = 0; i < N; ++i )
         {
            // find k that minimizes norm(x_i - mu_k)
            int min_k = -1;
            float min_dist;

            for ( int k = 0; k < K; ++k )
            {
               float dist = vectorDiffNorm(X[i], _means[k]);

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

         // update labels
         std::swap(y, y_next);

         // update means
         for ( int k = 0; k < K; ++k )
         {
            // compute mu_k = mean of all x_i in cluster k
            int n_k = 0;

            vectorInitZero(_means[k]);

            for ( int i = 0; i < N; ++i )
            {
               if ( y[i] == k )
               {
                  vectorAdd(_means[k], X[i]);
                  n_k++;
               }
            }

            vectorScale(_means[k], 1.0f / n_k);
         }
      }

      // save the run with the greatest log-likelihood
      float logL = computeLogLikelihood(X, N, y);

      if ( _logL < logL )
      {
         _logL = logL;
         std::swap(labels, y);
      }
   }

   return true;
}






float KMeans::computeLogLikelihood(const QVector<Vector2>& X, int N, const QVector<qint8>& y)
{
   // compute within-class scatter
   float S = 0;

   for ( int k = 0; k < _means.size(); ++k )
   {
      for ( int i = 0; i < N; ++i )
      {
         if ( y[i] != k )
         {
            continue;
         }

         float dist = vectorDiffNorm(X[i], _means[k]);

         S += dist * dist;
      }
   }

   return -S;
}
