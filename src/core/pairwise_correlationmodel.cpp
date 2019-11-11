#include "pairwise_correlationmodel.h"



using namespace Pairwise;



/*!
 * Compute the correlation of each cluster in a pairwise data array.
 *
 * @param expressions
 * @param index
 * @param K
 * @param labels
 * @param minSamples
 */
QVector<float> CorrelationModel::compute(
    const std::vector<float>& expressions,
    const Index& index,
    int K,
    const QVector<qint8>& labels,
    int minSamples)
{
    const float *x = &expressions[index.getX() * labels.size()];
    const float *y = &expressions[index.getY() * labels.size()];
    QVector<float> correlations(K);

    for ( qint8 k = 0; k < K; ++k )
    {
        correlations[k] = computeCluster(x, y, labels, k, minSamples);
    }

    return correlations;
}
