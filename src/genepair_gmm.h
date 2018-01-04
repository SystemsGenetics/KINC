#ifndef GENEPAIR_GMM_H
#define GENEPAIR_GMM_H
#include <ace/core/AceCore.h>

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

         void initialize(float pi, float *mu);
         void prepareCovariance();
         void calcLogMvNorm(const float *X, int N, int D, float *logP);

         float _pi;
         float *_mu;
         float *_sigma;

      private:
         float *_sigmaInv;
         float _normalizer;
      };

      int numClusters() const { return _K; }
      bool success() const { return _success; }
      float logLikelihood() const { return _logL; }
      const QVector<int>& labels() const { return _labels; }

      void fit(const float *X, int N, int D, int K, int maxIterations=100);

   private:
      void initialize(const float *X, int N, int D, int K);
      void kmeans(const float *X, int N, int D, const QVector<float *>& means);
      void calcLogMvNorm(const float *X, int N, int D, float *logProb);
      void calcLogLikelihoodAndGammaNK(const float *logpi, float *logProb, int N, float *logL);
      void calcLogGammaK(const float *loggamma, int N, float *logGamma);
      float calcLogGammaSum(const float *logpi, const float *logGamma);
      void performMStep(float *logpi, float *loggamma, float *logGamma, float logGammaSum, const float *X, int N, int D);
      QVector<int> calcLabels(float *loggamma, int N);

      int _K;
      QVector<Component> _components;
      bool _success;
      float _logL;
      QVector<int> _labels;
   };
}

#endif
