#ifndef EXTRACT_H
#define EXTRACT_H
#include <ace/core/AceCore.h>

#include "ccmatrix.h"
#include "correlationmatrix.h"
#include "expressionmatrix.h"



class Extract : public EAbstractAnalytic
{
   Q_OBJECT

public:
   enum Arguments
   {
      ExpressionData = 0
      ,ClusterData
      ,CorrelationData
      ,OutputFile
      ,GraphMLFile
      ,MinCorrelation
      ,MaxCorrelation
      ,Total
   };

   virtual int getArgumentCount() override final { return Total; }
   virtual ArgumentType getArgumentData(int argument) override final;
   virtual QVariant getArgumentData(int argument, Role role) override final;
   virtual void setArgument(int argument, EAbstractData* data) override final;
   virtual void setArgument(int argument, QFile* file) override final;
   virtual void setArgument(int argument, QVariant value) override final;
   quint32 getCapabilities() const override final { return Capabilities::Serial; }
   virtual bool initialize() override final;
   virtual void runSerial() override final;

private:
   ExpressionMatrix* _eMatrix {nullptr};
   CCMatrix* _ccMatrix {nullptr};
   CorrelationMatrix* _cMatrix {nullptr};
   QFile* _output {nullptr};
   QFile* _graphml {nullptr};
   float _minCorrelation {0.85};
   float _maxCorrelation {1.00};
};



#endif
