#include "genepair_kmeans.h"



using namespace GenePair;






void KMeans::fit(const QVector<Vector2>& X, int K)
{
   const int N = X.size();

   // initialize means randomly from X
   _means.resize(K);

   for ( int k = 0; k < K; ++k )
   {
      int i = rand() % N;
      _means[k] = X[i];
   }

   // iterate K means until convergence
   QVector<int> y;
   QVector<int> y_next(N);

   while ( true )
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
      y = y_next;

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

   // save outputs
   _logL = computeLogLikelihood(X);
   _labels = y;
}






float KMeans::computeLogLikelihood(const QVector<Vector2>& X)
{
   // compute within-class scatter
   float S = 0;

   for ( int k = 0; k < _means.size(); ++k )
   {
      for ( int i = 0; i < X.size(); ++i )
      {
         if ( _labels[i] != k )
         {
            continue;
         }

         float dist = vectorDiffNorm(X[i], _means[k]);

         S += dist * dist;
      }
   }

   return -S;
}
