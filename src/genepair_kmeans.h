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
      const QVector<cl_char>& labels() const { return _labels; }

      void fit(const QVector<Vector2>& X, int K, int numInits, int maxIterations);

   private:
      float computeLogLikelihood(const QVector<Vector2>& X, const QVector<cl_char>& y);

      QVector<Vector2> _means;
      float _logL;
      QVector<cl_char> _labels;
   };
}

#endif
