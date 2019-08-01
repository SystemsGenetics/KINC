#ifndef PAIRWISE_GMM_H
#define PAIRWISE_GMM_H
#include "pairwise_clusteringmodel.h"
#include "pairwise_linalg.h"

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
      virtual qint8 compute(
         const std::vector<float>& expressions,
         const Index& index,
         int numSamples,
         QVector<qint8>& labels,
         int minSamples,
         qint8 minClusters,
         qint8 maxClusters,
         Criterion criterion
      ) override final;
   private:
      void initializeComponents(const QVector<Vector2>& X, int N, int K);
      bool prepareComponents(int K);
      void initializeMeans(const QVector<Vector2>& X, int N, int K);
      float computeEStep(const QVector<Vector2>& X, int N, int K);
      void computeMStep(const QVector<Vector2>& X, int N, int K);
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
       * The array of mixture weights for each component.
       */
      float *_pi;
      /*!
       * The array of means for each component.
       */
      Vector2 *_mu;
      /*!
       * The array of covariance matrices for each component.
       */
      Matrix2x2 *_sigma;
      /*!
       * The array of precision matrices (inverse covariance) for each component.
       */
      Matrix2x2 *_sigmaInv;
      /*!
       * The array of normalizer terms for each component. The normalizer term
       * is pre-computed for the multivariate normal distribution function.
       */
      float *_normalizer;
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
