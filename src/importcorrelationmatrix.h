#ifndef IMPORTCORRELATIONMATRIX_H
#define IMPORTCORRELATIONMATRIX_H
#include <ace/core/AceCore.h>

#include "ccmatrix.h"
#include "correlationmatrix.h"



class ImportCorrelationMatrix : public EAbstractAnalytic
{
   Q_OBJECT

public:
   enum Arguments
   {
      InputFile = 0
      ,ClusterData
      ,CorrelationData
      ,GeneSize
      ,SampleSize
      ,CorrelationName
      ,Total
   };

   virtual int getArgumentCount() override final { return Total; }
   virtual ArgumentType getArgumentData(int argument) override final;
   virtual QVariant getArgumentData(int argument, Role role) override final;
   virtual void setArgument(int argument, QVariant value);
   virtual void setArgument(int argument, QFile* file) override final;
   virtual void setArgument(int argument, EAbstractData* data) override final;
   quint32 getCapabilities() const override final { return Capabilities::Serial; }
   virtual bool initialize() override final;
   virtual void runSerial() override final;

private:
   QFile* _input {nullptr};
   CCMatrix* _ccMatrix {nullptr};
   CorrelationMatrix* _cMatrix {nullptr};
   qint32 _geneSize {0};
   qint32 _sampleSize {0};
   QString _correlationName;
};



#endif
