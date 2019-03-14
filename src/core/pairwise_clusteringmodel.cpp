#include "pairwise_clusteringmodel.h"



using namespace Pairwise;






/*!
 * Construct an abstract pairwise clustering model.
 *
 * @param emx
 */
ClusteringModel::ClusteringModel(ExpressionMatrix* emx)
{
   // pre-allocate workspace
   _data.resize(emx->sampleSize());
   _labels.resize(emx->sampleSize());
}






/*!
 * Determine the number of clusters in a pairwise data array. Several sub-models,
 * each one having a different number of clusters, are fit to the data and the
 * sub-model with the best criterion value is selected.
 *
 * @param data
 * @param numSamples
 * @param labels
 * @param minSamples
 * @param minClusters
 * @param maxClusters
 * @param criterion
 */
qint8 ClusteringModel::compute(
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
      // extract clean samples from data array
      for ( int i = 0, j = 0; i < labels.size(); ++i )
      {
         if ( labels[i] >= 0 )
         {
            _data[j] = data[i];
            ++j;
         }
      }

      // determine the number of clusters
      float bestValue = INFINITY;

      for ( qint8 K = minClusters; K <= maxClusters; ++K )
      {
         // run each clustering sub-model
         bool success = fit(_data, numSamples, K, _labels);

         if ( !success )
         {
            continue;
         }

         // compute the criterion value of the sub-model
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

         // keep the sub-model with the lowest criterion value
         if ( value < bestValue )
         {
            bestK = K;
            bestValue = value;

            // save labels for clean samples
            for ( int i = 0, j = 0; i < labels.size(); ++i )
            {
               if ( labels[i] >= 0 )
               {
                  labels[i] = _labels[j];
                  ++j;
               }
            }
         }
      }
   }

   return bestK;
}
