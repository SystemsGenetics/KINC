#ifndef GENEPAIR_GMM_H
#define GENEPAIR_GNN_H
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

         void initialize(int D, float pi, float *mu);
         void prepareCovariance();
         void logMvNormDist(const float *X, int N, float *P);

         float _pi;
         float *_mu;
         float *_sigma;

      private:
         int _D;

         float *_sigmaL;
         float _normalizer;
      };

      int numClusters() const { return _K; }
      bool success() const { return _success; }
      float logLikelihood() const { return _logL; }
      const QVector<int>& labels() const { return _labels; }

      void fit(const float *X, int N, int D, int K, int maxIterations=100);

   private:
      void initialize(const float *X, int N, int D, int K);
      void kmeans(const float *X, int N, int D, float *M, int K);
      void calcLogMvNorm(const float *X, int N, float *logProb);
      void logLikelihoodAndGammaNK(const float *logpi, float *logProb, int N, float *logL);
      void calcLogGammaK(const float *loggamma, int N, float *logGamma);
      float calcLogGammaSum(const float *logpi, const float *logGamma);
      void performMStep(float *logpi, float *loggamma, float *logGamma, const float logGammaSum, const float *X, int N, int D, float *outerProduct, float *xm);
      QVector<int> calcLabels(float *loggamma, int N);

      int _K;
      QVector<Component> _components;
      bool _success;
      float _logL;
      QVector<int> _labels;
   };
}

#endif
