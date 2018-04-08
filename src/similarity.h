#ifndef SIMILARITY_H
#define SIMILARITY_H
#include <ace/core/AceCore.h>

#include "expressionmatrix.h"
#include "ccmatrix.h"
#include "correlationmatrix.h"
#include "genepair_clustering.h"
#include "genepair_correlation.h"



class Similarity : public EAbstractAnalytic
{
   Q_OBJECT

public:
   enum Arguments
   {
      InputData = 0
      ,ClusterData
      ,CorrelationData
      ,ClusteringArg
      ,CorrelationArg
      ,MinExpression
      ,MinSamples
      ,MinClusters
      ,MaxClusters
      ,CriterionArg
      ,RemovePreOutliers
      ,RemovePostOutliers
      ,MinCorrelation
      ,MaxCorrelation
      ,Total
   };

   enum class ClusteringMethod
   {
      GMM
      ,KMeans
   };

   enum class CorrelationMethod
   {
      Pearson
      ,Spearman
   };

   virtual int getArgumentCount() override final { return Total; }
   virtual ArgumentType getArgumentData(int argument) override final;
   virtual QVariant getArgumentData(int argument, Role role) override final;
   virtual void setArgument(int argument, QVariant value) override final;
   virtual void setArgument(int argument, EAbstractData* data) override final;
   quint32 getCapabilities() const override final
      { return Capabilities::Serial; }
   virtual bool initialize() override final;
   virtual void runSerial() override final;
   virtual void finish() override final {}

private:
   static const char* GMM;
   static const char* KMeans;
   static const char* Pearson;
   static const char* Spearman;
   static const char* BIC;
   static const char* ICL;

   void savePair(GenePair::Vector vector, qint8 K, const qint8 *labels, int N, const float *correlations);

   ExpressionMatrix* _input {nullptr};
   CCMatrix* _clusMatrix {nullptr};
   CorrelationMatrix* _corrMatrix {nullptr};
   ClusteringMethod _clusMethod {ClusteringMethod::GMM};
   CorrelationMethod _corrMethod {CorrelationMethod::Pearson};
   GenePair::Clustering* _clusModel {nullptr};
   GenePair::Correlation* _corrModel {nullptr};
   int _minSamples {30};
   float _minExpression {-INFINITY};
   qint8 _minClusters {1};
   qint8 _maxClusters {5};
   GenePair::Criterion _criterion {GenePair::Criterion::BIC};
   bool _removePreOutliers {false};
   bool _removePostOutliers {false};
   float _minCorrelation {0.5};
   float _maxCorrelation {1.0};
};



#endif
