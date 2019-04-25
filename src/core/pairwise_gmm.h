#ifndef PAIRWISE_GMM_H
#define PAIRWISE_GMM_H
#include "pairwise_clusteringmodel.h"

namespace Pairwise
{
   /*!
    * This class implements the Gaussian mixture model. The number of clusters is
    * determined by creating several sub-models, each with a different assumption
    * of the number of clusters, and selecting the sub-model which best fits the
    * data according to a criterion.
    */
   class GMM : public ClusteringModel
   {
   public:
      GMM(ExpressionMatrix* emx, qint8 maxClusters);
      ~GMM();
   public:
      class Component
      {
      public:
         Component() = default;
         void initialize(float pi, const Vector2& mu);
         bool prepare();
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
   public:
      virtual qint8 compute(
         const QVector<Vector2>& data,
         int numSamples,
         QVector<qint8>& labels,
         int minSamples,
         qint8 minClusters,
         qint8 maxClusters,
         Criterion criterion
      ) override final;
   private:
      void initializeMeans(const QVector<Vector2>& X, int N);
      float computeEStep(const QVector<Vector2>& X, int N);
      void computeMStep(const QVector<Vector2>& X, int N);
      void computeLabels(const float *gamma, int N, int K, QVector<qint8>& labels);
      float computeEntropy(const float *gamma, int N, const QVector<qint8>& labels);
      bool fit(const QVector<Vector2>& X, int N, int K, QVector<qint8>& labels);
      float computeAIC(int K, int D, float logL);
      float computeBIC(int K, int D, float logL, int N);
      float computeICL(int K, int D, float logL, int N, float E);
   private:
      /*!
       * Workspace for clustering data.
       */
      QVector<Vector2> _data;
      /*!
       * Workspace for the cluster labels.
       */
      QVector<qint8> _labels;
      /*!
       * The list of mixture components, which define the mean and covariance
       * of each cluster in the mixture model.
       */
      QVector<Component> _components;
      /*!
       * The array of posterior probabilities used by the EM algorithm.
       */
      float *_gamma;
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
