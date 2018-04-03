#ifndef GENEPAIR_KMEANS_H
#define GENEPAIR_KMEANS_H
#include <ace/core/AceCore.h>

#include "genepair_clustering.h"

namespace GenePair
{
   class KMeans : public Clustering
   {
   public:
      KMeans() = default;

      bool fit(const QVector<Vector2>& X, int K, QVector<qint8>& labels);

      float logLikelihood() const { return _logL; }
      float entropy() const { return 0; }

   private:
      float computeLogLikelihood(const QVector<Vector2>& X, const QVector<qint8>& y);

      QVector<Vector2> _means;
      float _logL;
   };
}

#endif
