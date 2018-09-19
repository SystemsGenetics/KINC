#include "pairwise_clustering.h"



using namespace Pairwise;






Clustering::Clustering(ExpressionMatrix* input)
{
   // pre-allocate workspace
   _workLabels.resize(input->sampleSize());
}






/*!
 * Compute clusters for a given dataset. Several clustering models, each one
 * having a different number of clusters, are fit to the data and the model
 * with the best criterion value is selected.
 *
 * Note that the dataset contains only those samples which were not removed
 * by pre-processing, while the labels contains all samples from the original
 * expression matrix.
 *
 * @param data
 * @param numSamples
 * @param labels
 * @param minSamples
 * @param minSamples
 * @param minClusters
 * @param maxClusters
 * @param criterion
 */
qint8 Clustering::compute(
   const QVector<Vector2>& data,
   int numSamples,
   QVector<qint8>& labels,
   int minSamples,
   qint8 minClusters,
   qint8 maxClusters,
   Criterion criterion)
{
   // perform clustering only if there are enough samples
   qint8 bestK = 0;

   if ( numSamples >= minSamples )
   {
      float bestValue = INFINITY;

      for ( qint8 K = minClusters; K <= maxClusters; ++K )
      {
         // run each clustering model
         bool success = fit(data, numSamples, K, _workLabels);

         if ( !success )
         {
            continue;
         }

         // evaluate model
         float value = INFINITY;

         switch (criterion)
         {
         case Criterion::AIC:
            value = computeAIC(K, 2, logLikelihood());
            break;
         case Criterion::BIC:
            value = computeBIC(K, 2, logLikelihood(), numSamples);
            break;
         case Criterion::ICL:
            value = computeICL(K, 2, logLikelihood(), numSamples, entropy());
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

   return bestK;
}






float Clustering::computeAIC(int K, int D, float logL)
{
   int p = K * (1 + D + D * D);

   return 2 * p - 2 * logL;
}






float Clustering::computeBIC(int K, int D, float logL, int N)
{
   int p = K * (1 + D + D * D);

   return log(N) * p - 2 * logL;
}






float Clustering::computeICL(int K, int D, float logL, int N, float E)
{
   int p = K * (1 + D + D * D);

   return log(N) * p - 2 * logL - 2 * E;
}
