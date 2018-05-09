#include "genepair_correlation.h"



using namespace Pairwise;






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
