#include "pairwise_correlationmodel.h"



using namespace Pairwise;






/*!
 * Compute the correlation of each cluster in a pairwise data array.
 *
 * @param data
 * @param K
 * @param labels
 * @param minSamples
 */
QVector<float> CorrelationModel::compute(
   const QVector<Vector2>& data,
   int K,
   const QVector<qint8>& labels,
   int minSamples)
{
   QVector<float> correlations(K);

   for ( qint8 k = 0; k < K; ++k )
   {
      correlations[k] = computeCluster(data, labels, k, minSamples);
   }

   return correlations;
}
