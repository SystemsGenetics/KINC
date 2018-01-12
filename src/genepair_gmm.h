#ifndef GENEPAIR_GMM_H
#define GENEPAIR_GMM_H
#include <ace/core/AceCore.h>

#include "genepair_linalg.h"

namespace GenePair
{
   class GMM
   {
   public:
      GMM() = default;

      class Component
      {
      public:
         Component() = default;

         void initialize(float pi, const Vector2& mu);
         void prepareCovariance();
         void calcLogMvNorm(const QVector<Vector2>& X, float *logP);

         float _pi;
         Vector2 _mu;
         Matrix2x2 _sigma;

      private:
         Matrix2x2 _sigmaInv;
         float _normalizer;
      };

      int numClusters() const { return _components.size(); }
      bool success() const { return _success; }
      float logLikelihood() const { return _logL; }
      float entropy() const { return _entropy; }
      const QVector<int>& labels() const { return _labels; }

      void fit(const QVector<Vector2>& X, int K);

   private:
      void kmeans(const QVector<Vector2>& X);
      void calcLogMvNorm(const QVector<Vector2>& X, float *loggamma);
      void calcLogLikelihoodAndGammaNK(const float *logpi, int K, float *loggamma, int N, float *logL);
      void calcLogGammaK(const float *loggamma, int N, int K, float *logGamma);
      float calcLogGammaSum(const float *logpi, int K, const float *logGamma);
      void performMStep(float *logpi, int K, float *loggamma, float *logGamma, float logGammaSum, const QVector<Vector2>& X);
      QVector<int> calcLabels(float *loggamma, int N, int K);
      float calcEntropy(float *loggamma, int N);

      QVector<Component> _components;
      bool _success;
      float _logL;
      float _entropy;
      QVector<int> _labels;
   };
}

#endif
