#ifndef GENEPAIR_KMEANS_H
#define GENEPAIR_KMEANS_H
#include <ace/core/AceCore.h>

#include "genepair_linalg.h"

namespace GenePair
{
   class KMeans
   {
   public:
      KMeans() = default;

      int numClusters() const { return _means.size(); }
      float logLikelihood() const { return _logL; }
      const QVector<int>& labels() const { return _labels; }

      void fit(const QVector<Vector2>& X, int K);

   private:
      void initialize(const QVector<Vector2>& X, int K);
      float computeLogLikelihood(const QVector<Vector2>& X);

      QVector<Vector2> _means;
      float _logL;
      QVector<int> _labels;
   };
}

#endif
