#ifndef PAIRWISE_CLUSTERINGMODEL_H
#define PAIRWISE_CLUSTERINGMODEL_H
#include <ace/core/core.h>

#include "ccmatrix.h"
#include "expressionmatrix.h"
#include "pairwise_index.h"

namespace Pairwise
{
    /*!
     * Defines the criterion types used by the abstract clustering model.
     */
    enum class Criterion
    {
        /*!
         * Akaike information criterion
         */
        AIC
        /*!
         * Bayesian information criterion
         */
        ,BIC
        /*!
         * Integrated completed likelihood
         */
        ,ICL
    };

    /*!
     * This class implements the abstract pairwise clustering model, which takes
     * a pairwise data array and determines the number of clusters, as well as the
     * cluster label for each sample in the data array. The underlying clustering
     * model must be implemented by the inheriting class.
     */
    class ClusteringModel
    {
    public:
        ~ClusteringModel() = default;
    public:
        virtual qint8 compute(
            const std::vector<float>& expressions,
            const Index& index,
            int numSamples,
            QVector<qint8>& labels,
            int minSamples,
            qint8 minClusters,
            qint8 maxClusters,
            Criterion criterion
        ) = 0;
    };
}

#endif
