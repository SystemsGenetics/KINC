#ifndef EXPORTEXPRESSIONMATRIX_H
#define EXPORTEXPRESSIONMATRIX_H
#include <ace/core/AceCore.h>

#include "expressionmatrix.h"



class ExportExpressionMatrix : public EAbstractAnalytic
{
   Q_OBJECT
public:
   enum Arguments
   {
      InputData = 0
      ,OutputFile
      ,NoSampleToken
      ,Total
   };
   virtual int getArgumentCount() override final { return Total; }
   virtual ArgumentType getArgumentData(int argument) override final;
   virtual QVariant getArgumentData(int argument, Role role) override final;
   virtual void setArgument(int argument, QVariant value) override final;
   virtual void setArgument(int argument, EAbstractData* data) override final;
   virtual void setArgument(int argument, QFile* file) override final;
   quint32 getCapabilities() const override final { return Capabilities::Serial; }
   virtual bool initialize() override final;
   virtual void runSerial() override final;
private:
   ExpressionMatrix* _input {nullptr};
   QFile* _output {nullptr};
   QString _noSampleToken;
};



#endif
