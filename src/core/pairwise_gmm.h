#ifndef PAIRWISE_GMM_H
#define PAIRWISE_GMM_H
#include "pairwise_clusteringmodel.h"

namespace Pairwise
{
   /*!
    * This class implements the Gaussian mixture model.
    */
   class GMM : public ClusteringModel
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
         void computeLogProbNorm(const QVector<Vector2>& X, int N, float *logP);
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
      void initializeMeans(const QVector<Vector2>& X, int N);
      float computeEStep(const QVector<Vector2>& X, int N, const float *logpi, float *loggamma, float *logGamma, float *logGammaSum);
      void computeMStep(const QVector<Vector2>& X, int N, float *logpi, float *loggamma, float *logGamma, float logGammaSum);
      void computeLabels(const float *gamma, int N, int K, QVector<qint8>& labels);
      float computeEntropy(const float *gamma, int N, const QVector<qint8>& labels);
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
