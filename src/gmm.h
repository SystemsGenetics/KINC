#ifndef GMM_H
#define GMM_H
#include <ace/core/AceCore.h>

#include "expressionmatrix.h"
#include "ccmatrix.h"
#include "genepair_gmm.h"
#include "genepair_vector.h"



class GMM : public EAbstractAnalytic
{
   Q_OBJECT

public:
   enum Arguments
   {
      InputData = 0
      ,OutputData
      ,MinSamples
      ,MinClusters
      ,MaxClusters
      ,Total
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
   float computeBIC(const GenePair::GMM& gmm, int N, int D);
   void computeModel(const QVector<GenePair::Vector2>& X, int& bestK, QVector<int>& bestLabels);
   void savePair(int K, const QVector<int>& labels, const GenePair::Vector& vector);

   ExpressionMatrix* _input {nullptr};
   CCMatrix* _output {nullptr};
   int _minSamples {30};
   int _minClusters {1};
   int _maxClusters {5};
};



#endif
