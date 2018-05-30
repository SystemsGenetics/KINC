#ifndef SIMILARITY_H
#define SIMILARITY_H
#include <ace/core/core.h>

#include "ccmatrix.h"
#include "correlationmatrix.h"
#include "expressionmatrix.h"
#include "pairwise_clustering.h"
#include "pairwise_correlation.h"
#include "pairwise_gmm.h"
#include "pairwise_pearson.h"



class Similarity : public EAbstractAnalytic
{
   Q_OBJECT
public:
   struct Pair
   {
      qint8 K;
      QVector<qint8> labels;
      QVector<float> correlations;
   };

   class Input;
   class WorkBlock;
   class ResultBlock;
   class Serial;
   class OpenCL;
   virtual int size() const override final;
   virtual std::unique_ptr<EAbstractAnalytic::Block> makeWork(int index) const override final;
   virtual std::unique_ptr<EAbstractAnalytic::Block> makeWork() const override final;
   virtual std::unique_ptr<EAbstractAnalytic::Block> makeResult() const override final;
   virtual void process(const EAbstractAnalytic::Block* result) override final;
   virtual EAbstractAnalytic::Input* makeInput() override final;
   virtual EAbstractAnalytic::Serial* makeSerial() override final;
   virtual EAbstractAnalytic::OpenCL* makeOpenCL() override final;
   virtual void initialize() override final;

private:
   enum class ClusteringMethod
   {
      None
      ,GMM
      ,KMeans
   };

   enum class CorrelationMethod
   {
      Pearson
      ,Spearman
   };

   ExpressionMatrix* _input {nullptr};
   CCMatrix* _ccm {nullptr};
   CorrelationMatrix* _cmx {nullptr};
   ClusteringMethod _clusMethod {ClusteringMethod::None};
   CorrelationMethod _corrMethod {CorrelationMethod::Pearson};
   Pairwise::Clustering* _clusModel {nullptr};
   Pairwise::Correlation* _corrModel {new Pairwise::Pearson()};
   int _minSamples {30};
   float _minExpression {-std::numeric_limits<float>::infinity()};
   qint8 _minClusters {1};
   qint8 _maxClusters {5};
   Pairwise::Criterion _criterion {Pairwise::Criterion::ICL};
   bool _removePreOutliers {false};
   bool _removePostOutliers {false};
   float _minCorrelation {0.5};
   float _maxCorrelation {1.0};
   int _kernelSize {4096};
};



#endif
