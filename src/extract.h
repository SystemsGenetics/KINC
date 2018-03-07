#ifndef EXTRACT_H
#define EXTRACT_H
#include <ace/core/AceCore.h>

#include "ccmatrix.h"
#include "correlationmatrix.h"



class Extract : public EAbstractAnalytic
{
   Q_OBJECT

public:
   enum Arguments
   {
      ClusterData = 0
      ,CorrelationData
      ,OutputFile
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
   CCMatrix* _ccMatrix {nullptr};
   CorrelationMatrix* _cMatrix {nullptr};
   QFile* _output {nullptr};
   float _minCorrelation {0.85};
   float _maxCorrelation {1.00};
};



#endif
