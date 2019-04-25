#ifndef SIMILARITY_H
#define SIMILARITY_H
#include <ace/core/core.h>

#include "ccmatrix.h"
#include "correlationmatrix.h"
#include "expressionmatrix.h"
#include "pairwise_clusteringmodel.h"



/*!
 * This class implements the similarity analytic. This analytic takes an
 * expression matrix and computes a similarity matrix, where each element is
 * a similarity measure of two genes in the expression matrix. The similarity
 * is computed using a correlation measure. The similarity matrix can also have
 * multiple modes within a pair; these modes can be optionally computed using a
 * clustering method. This analytic produces two data objects: a correlation
 * matrix containing the pairwise correlations, and a cluster matrix containing
 * sample masks of the pairwise clusters. Sample masks for unimodal pairs are not
 * saved to the cluster matrix. If clustering is not used, an empty cluster matrix
 * is created. This analytic can also perform pairwise outlier removal before and
 * after clustering, if clustering is used.
 *
 * This analytic can use MPI and it has both CPU and GPU implementations, as the
 * pairwise clustering significantly increases the amount of computations required
 * for a large expression matrix.
 */
class Similarity : public EAbstractAnalytic
{
   Q_OBJECT
public:
   /*!
    * Defines the pair structure used to send results in result blocks.
    */
   struct Pair
   {
      /*!
       * The number of clusters in a pair.
       */
      qint8 K;
      /*!
       * The cluster labels for a pair.
       */
      QVector<qint8> labels;
      /*!
       * The correlation for each cluster in a pair.
       */
      QVector<float> correlations;
   };
   class Input;
   class WorkBlock;
   class ResultBlock;
   class Serial;
   class OpenCL;
   class CUDA;
public:
   static int nextPower2(int n);
   static qint64 totalPairs(const ExpressionMatrix* emx);
public:
   virtual int size() const override final;
   virtual std::unique_ptr<EAbstractAnalyticBlock> makeWork(int index) const override final;
   virtual std::unique_ptr<EAbstractAnalyticBlock> makeWork() const override final;
   virtual std::unique_ptr<EAbstractAnalyticBlock> makeResult() const override final;
   virtual void process(const EAbstractAnalyticBlock* result) override final;
   virtual EAbstractAnalyticInput* makeInput() override final;
   virtual EAbstractAnalyticSerial* makeSerial() override final;
   virtual EAbstractAnalyticOpenCL* makeOpenCL() override final;
   virtual EAbstractAnalyticCUDA* makeCUDA() override final;
   virtual void initialize() override final;
   virtual void initializeOutputs() override final;
private:
   /*!
    * Defines the clustering methods this analytic supports.
    */
   enum class ClusteringMethod
   {
      /*!
       * No clustering
       */
      None
      /*!
       * Gaussian mixture models
       */
      ,GMM
   };
   /*!
    * Defines the correlation methods this analytic supports.
    */
   enum class CorrelationMethod
   {
      /*!
       * Pearson correlation
       */
      Pearson
      /*!
       * Spearman rank correlation
       */
      ,Spearman
   };
private:
   /*!
    * Pointer to the input expression matrix.
    */
   ExpressionMatrix* _input {nullptr};
   /*!
    * Pointer to the output cluster matrix.
    */
   CCMatrix* _ccm {nullptr};
   /*!
    * Pointer to the output correlation matrix.
    */
   CorrelationMatrix* _cmx {nullptr};
   /*!
    * The clustering method to use.
    */
   ClusteringMethod _clusMethod {ClusteringMethod::None};
   /*!
    * The correlation method to use.
    */
   CorrelationMethod _corrMethod {CorrelationMethod::Pearson};
   /*!
    * The name of the correlation method.
    */
   QString _corrName;
   /*!
    * The minimum number of clean samples required to consider a pair.
    */
   int _minSamples {30};
   /*!
    * The minimum expression value required to include a sample.
    */
   float _minExpression {-std::numeric_limits<float>::infinity()};
   /*!
    * The minimum number of clusters to use in the clustering model.
    */
   qint8 _minClusters {1};
   /*!
    * The maximum number of clusters to use in the clustering model.
    */
   qint8 _maxClusters {5};
   /*!
    * The model selection criterion to use in the clustering model.
    */
   Pairwise::Criterion _criterion {Pairwise::Criterion::ICL};
   /*!
    * Whether to remove outliers before clustering.
    */
   bool _removePreOutliers {false};
   /*!
    * Whether to remove outliers after clustering.
    */
   bool _removePostOutliers {false};
   /*!
    * The minimum (absolute) correlation threshold to save a correlation.
    */
   float _minCorrelation {0.5};
   /*!
    * The maximum (absolute) correlation threshold to save a correlation.
    */
   float _maxCorrelation {1.0};
   /*!
    * The number of pairs to process in each work block.
    */
   int _workBlockSize {0};
   /*!
    * The global work size for each OpenCL worker.
    */
   int _globalWorkSize {4096};
   /*!
    * The local work size for each OpenCL worker.
    */
   int _localWorkSize {32};
};



#endif
