#ifndef PAIRWISE_CLUSTERING_H
#define PAIRWISE_CLUSTERING_H
#include <ace/core/core.h>

#include "ccmatrix.h"
#include "expressionmatrix.h"
#include "pairwise_linalg.h"
#include "pairwise_index.h"

namespace Pairwise
{
   enum class Criterion
   {
      AIC
      ,BIC
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
         Criterion criterion
      );

   protected:
      virtual bool fit(const QVector<Vector2>& X, int N, int K, QVector<qint8>& labels) = 0;
      virtual float logLikelihood() const = 0;
      virtual float entropy() const = 0;

   private:
      float computeAIC(int K, int D, float logL);
      float computeBIC(int K, int D, float logL, int N);
      float computeICL(int K, int D, float logL, int N, float E);

      QVector<qint8> _workLabels;
   };
}

#endif
