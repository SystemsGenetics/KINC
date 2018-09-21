#ifndef PAIRWISE_GMM_H
#define PAIRWISE_GMM_H
#include "pairwise_clustering.h"

namespace Pairwise
{
   /*!
    * This class implements the Gaussian mixture model.
    */
   class GMM : public Clustering
   {
   public:
      GMM(ExpressionMatrix* emx);
   public:
      class Component
      {
      public:
         Component() = default;
         void initialize(float pi, const Vector2& mu);
         void prepare();
         void calcLogMvNorm(const QVector<Vector2>& X, int N, float *logP);
      public:
         /*!
          * The mixture weight.
          */
         float _pi;
         /*!
          * The mean.
          */
         Vector2 _mu;
         /*!
          * The covariance matrix.
          */
         Matrix2x2 _sigma;
      private:
         /*!
          * The precision matrix, or inverse of the covariance matrix.
          */
         Matrix2x2 _sigmaInv;
         /*!
          * A normalization term which is pre-computed for the multivariate
          * normal distribution function.
          */
         float _normalizer;
      };
   protected:
      bool fit(const QVector<Vector2>& X, int N, int K, QVector<qint8>& labels);
      float logLikelihood() const { return _logL; }
      float entropy() const { return _entropy; }
      float computeAIC(int K, int D, float logL);
      float computeBIC(int K, int D, float logL, int N);
      float computeICL(int K, int D, float logL, int N, float E);
   private:
      void kmeans(const QVector<Vector2>& X, int N);
      void calcLogMvNorm(const QVector<Vector2>& X, int N, float *loggamma);
      void calcLogLikelihoodAndGammaNK(const float *logpi, int K, float *loggamma, int N, float *logL);
      void calcLogGammaK(const float *loggamma, int N, int K, float *logGamma);
      float calcLogGammaSum(const float *logpi, int K, const float *logGamma);
      void performMStep(float *logpi, int K, float *loggamma, float *logGamma, float logGammaSum, const QVector<Vector2>& X, int N);
      void calcLabels(float *loggamma, int N, int K, QVector<qint8>& labels);
      float calcEntropy(float *loggamma, int N, const QVector<qint8>& labels);
      /*!
       * The list of mixture components, which define the mean and covariance
       * of each cluster in the mixture model.
       */
      QVector<Component> _components;
      /*!
       * The log-likelihood of the mixture model.
       */
      float _logL;
      /*!
       * The entropy of the mixture model.
       */
      float _entropy;
   };
}

#endif
