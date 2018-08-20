#include "pairwise_correlation.h"



using namespace Pairwise;






/*!
 * Compute cluster-wise correlations for a given dataset and labels.
 *
 * Note that the dataset contains only those samples which were not removed
 * by pre-processing, while the labels contains all samples from the original
 * expression matrix.
 *
 * @param data
 * @param K
 * @param labels
 * @param minSamples
 */
QVector<float> Correlation::compute(
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
