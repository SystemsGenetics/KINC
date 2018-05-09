#ifndef PAIRWISE_CLUSTERING_H
#define PAIRWISE_CLUSTERING_H
#include <ace/core/AceCore.h>

#include "ccmatrix.h"
#include "expressionmatrix.h"
#include "genepair_linalg.h"
#include "genepair_index.h"

namespace Pairwise
{
   enum class Criterion
   {
      BIC
      ,ICL
   };

   class Clustering
   {
   public:
      void initialize(ExpressionMatrix* input);
      qint8 compute(
         const QVector<Vector2>& X,
         int numSamples,
         QVector<qint8>& labels,
         int minSamples,
         qint8 minClusters,
         qint8 maxClusters,
         Criterion criterion,
         bool removePreOutliers,
         bool removePostOutliers
      );

   protected:
      virtual bool fit(const QVector<Vector2>& X, int N, int K, QVector<qint8>& labels) = 0;
      virtual float logLikelihood() const = 0;
      virtual float entropy() const = 0;

   private:
      void markOutliers(const QVector<Vector2>& X, int N, int j, QVector<qint8>& labels, qint8 cluster, qint8 marker);
      float computeBIC(int K, float logL, int N, int D);
      float computeICL(int K, float logL, int N, int D, float E);

      QVector<qint8> _workLabels;
   };
}

#endif
