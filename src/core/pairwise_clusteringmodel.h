#ifndef PAIRWISE_CLUSTERINGMODEL_H
#define PAIRWISE_CLUSTERINGMODEL_H
#include <ace/core/core.h>

#include "ccmatrix.h"
#include "expressionmatrix.h"
#include "pairwise_linalg.h"
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
    * cluster label for each sample in the data array. The number of clusters is
    * determined by creating several sub-models, each with a different assumption
    * of the number of clusters, and selecting the sub-model which best fits the
    * data according to a criterion. The clustering sub-model must be implemented
    * by the inheriting class.
    */
   class ClusteringModel
   {
   public:
      ClusteringModel(ExpressionMatrix* emx);
      ~ClusteringModel() = default;
      qint8 compute(
         const QVector<Vector2>& data,
         int numSamples,
         QVector<qint8>& labels,
         int minSamples,
         qint8 minClusters,
         qint8 maxClusters,
         Criterion criterion
      );
   protected:
      virtual bool fit(const QVector<Vector2>& X, int N, int K, QVector<qint8>& labels) = 0;
      virtual float logLikelihood() const = 0;
      virtual float entropy() const = 0;
      virtual float computeAIC(int K, int D, float logL) = 0;
      virtual float computeBIC(int K, int D, float logL, int N) = 0;
      virtual float computeICL(int K, int D, float logL, int N, float E) = 0;
   private:
      /*!
       * Workspace for clustering data.
       */
      QVector<Vector2> _data;
      /*!
       * Workspace for the cluster labels.
       */
      QVector<qint8> _labels;
   };
}

#endif
