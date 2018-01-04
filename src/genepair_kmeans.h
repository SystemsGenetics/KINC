#ifndef GENEPAIR_KMEANS_H
#define GENEPAIR_KMEANS_H
#include <ace/core/AceCore.h>

namespace GenePair
{
   class KMeans
   {
   public:
      KMeans() = default;
      ~KMeans();

      int numClusters() const { return _K; }
      float logLikelihood() const { return _logL; }
      const QVector<int>& labels() const { return _labels; }

      void fit(const float *X, int N, int D, int K);

   private:
      void initialize(const float *X, int N, int D, int K);
      float computeLogLikelihood(const float *X, int N, int D);

      int _K;
      QVector<float *> _means;
      float _logL;
      QVector<int> _labels;
   };
}

#endif
