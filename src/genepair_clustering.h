#ifndef GENEPAIR_CLUSTERING_H
#define GENEPAIR_CLUSTERING_H
#include <ace/core/AceCore.h>

#include "expressionmatrix.h"
#include "genepair_linalg.h"
#include "genepair_vector.h"

namespace GenePair
{
   enum class Criterion
   {
      BIC
      ,ICL
   };

   class Clustering
   {
   public:
      void compute(
         ExpressionMatrix* input,
         Vector vector,
         int minSamples,
         int minExpression,
         qint8 minClusters,
         qint8 maxClusters,
         Criterion criterion,
         bool removePreOutliers,
         bool removePostOutliers
      );

      qint8 clusterSize() const { return _bestK; };
      const QVector<qint8>& labels() const { return _bestLabels; };

   protected:
      virtual bool fit(const QVector<Vector2>& X, int K, QVector<qint8>& labels) = 0;
      virtual float logLikelihood() const = 0;
      virtual float entropy() const = 0;

   private:
      void fetchData(ExpressionMatrix* input, Vector vector, int minExpression);
      void markOutliers(const QVector<Vector2>& X, int j, QVector<qint8>& labels, qint8 cluster, qint8 marker);
      float computeBIC(int K, float logL, int N, int D);
      float computeICL(int K, float logL, int N, int D, float E);

      QVector<Vector2> _X;
      QVector<qint8> _labels;
      qint8 _bestK;
      QVector<qint8> _bestLabels;
   };
}

#endif
