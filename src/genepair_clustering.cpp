#include "genepair_clustering.h"



using namespace Pairwise;






void Clustering::initialize(ExpressionMatrix* input)
{
   // pre-allocate workspace
   _workLabels.resize(input->getSampleSize());
}






qint8 Clustering::compute(
   const QVector<Vector2>& X,
   int numSamples,
   QVector<qint8>& labels,
   int minSamples,
   qint8 minClusters,
   qint8 maxClusters,
   Criterion criterion,
   bool removePreOutliers,
   bool removePostOutliers)
{
   // remove pre-clustering outliers
   if ( removePreOutliers )
   {
      markOutliers(X, numSamples, 0, labels, 0, -7);
      markOutliers(X, numSamples, 1, labels, 0, -7);
   }

   // perform clustering only if there are enough samples
   qint8 bestK = 0;

   if ( numSamples >= minSamples )
   {
      float bestValue = INFINITY;

      for ( qint8 K = minClusters; K <= maxClusters; ++K )
      {
         // run each clustering model
         bool success = fit(X, numSamples, K, _workLabels);

         if ( !success )
         {
            continue;
         }

         // evaluate model
         float value = INFINITY;

         switch (criterion)
         {
         case Criterion::BIC:
            value = computeBIC(K, logLikelihood(), numSamples, 2);
            break;
         case Criterion::ICL:
            value = computeICL(K, logLikelihood(), numSamples, 2, entropy());
            break;
         }

         // save the best model
         if ( value < bestValue )
         {
            bestK = K;
            bestValue = value;

            for ( int i = 0, j = 0; i < numSamples; ++i )
            {
               if ( labels[i] >= 0 )
               {
                  labels[i] = _workLabels[j];
                  ++j;
               }
            }
         }
      }
   }

   if ( bestK > 1 )
   {
      // remove post-clustering outliers
      if ( removePostOutliers )
      {
         for ( qint8 k = 0; k < bestK; ++k )
         {
            markOutliers(X, numSamples, 0, labels, k, -8);
            markOutliers(X, numSamples, 1, labels, k, -8);
         }
      }
   }

   return bestK;
}






void Clustering::markOutliers(const QVector<Vector2>& X, int N, int j, QVector<qint8>& labels, qint8 cluster, qint8 marker)
{
   // compute x_sorted = X[:, j], filtered and sorted
   QVector<float> x_sorted;
   x_sorted.reserve(N);

   for ( int i = 0; i < N; i++ )
   {
      if ( labels[i] == cluster || labels[i] == marker )
      {
         x_sorted.append(X[i].s[j]);
      }
   }

   if ( x_sorted.size() == 0 )
   {
      return;
   }

   std::sort(x_sorted.begin(), x_sorted.end());

   // compute quartiles, interquartile range, upper and lower bounds
   const int n = x_sorted.size();

   float Q1 = x_sorted[n * 1 / 4];
   float Q3 = x_sorted[n * 3 / 4];

   float T_min = Q1 - 1.5f * (Q3 - Q1);
   float T_max = Q3 + 1.5f * (Q3 - Q1);

   // mark outliers
   for ( int i = 0; i < N; ++i )
   {
      if ( labels[i] == cluster && (X[i].s[j] < T_min || T_max < X[i].s[j]) )
      {
         labels[i] = marker;
      }
   }
}






float Clustering::computeBIC(int K, float logL, int N, int D)
{
   int p = K * (1 + D + D * D);

   return log(N) * p - 2 * logL;
}






float Clustering::computeICL(int K, float logL, int N, int D, float E)
{
   int p = K * (1 + D + D * D);

   return log(N) * p - 2 * logL - 2 * E;
}
