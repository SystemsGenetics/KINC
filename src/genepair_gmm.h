#ifndef GENEPAIR_GMM_H
#define GENEPAIR_GMM_H
#include <ace/core/AceCore.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace GenePair
{
   class GMM
   {
   public:
      GMM();

      class Component
      {
      public:
         Component();
         ~Component();

         void initialize(float pi, gsl_vector_float *mu);
         void prepareCovariance();
         void logMvNormDist(const gsl_matrix_float *X, float *logProbK);

         float _pi;
         gsl_vector_float *_mu;
         gsl_matrix_float *_sigma;

      private:
         gsl_matrix_float *_sigmaInv;
         float _normalizer;
      };

      int numClusters() const { return _K; }
      bool success() const { return _success; }
      float logLikelihood() const { return _logL; }
      const QVector<int>& labels() const { return _labels; }

      void fit(const gsl_matrix_float *X, int K, int maxIterations=100);

   private:
      void initialize(const gsl_matrix_float *X, int K);
      void kmeans(const gsl_matrix_float *X, const QVector<gsl_vector_float *>& means);
      void calcLogMvNorm(const gsl_matrix_float *X, float *logProb);
      void logLikelihoodAndGammaNK(const float *logpi, float *logProb, int N, float *logL);
      void calcLogGammaK(const float *loggamma, int N, float *logGamma);
      float calcLogGammaSum(const float *logpi, const float *logGamma);
      void performMStep(float *logpi, float *loggamma, float *logGamma, float logGammaSum, const gsl_matrix_float *X);
      QVector<int> calcLabels(float *loggamma, int N);

      int _K;
      QVector<Component> _components;
      bool _success;
      float _logL;
      QVector<int> _labels;
   };
}

#endif
