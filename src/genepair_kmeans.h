#ifndef GENEPAIR_KMEANS_H
#define GENEPAIR_KMEANS_H
#include <ace/core/AceCore.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

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

      void fit(const gsl_matrix_float *X, int K);

   private:
      void initialize(const gsl_matrix_float *X, int K);
      float computeLogLikelihood(const gsl_matrix_float *X);

      int _K;
      QVector<gsl_vector_float *> _means;
      float _logL;
      QVector<int> _labels;
   };
}

#endif
