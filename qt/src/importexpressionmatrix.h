#ifndef IMPORTEXPRESSIONMATRIX_H
#define IMPORTEXPRESSIONMATRIX_H
#include <ace/core/AceCore.h>

#include "expressionmatrix.h"



class ImportExpressionMatrix : public EAbstractAnalytic
{
public:
   enum Arguments
   {
      InputFile = 0
      ,OutputData
      ,NoSampleToken
      ,SampleSize
      ,TransformArg
      ,Total
   };
   virtual int getArgumentCount() override final { return Total; }
   virtual ArgumentType getArgumentData(int argument) override final;
   virtual QVariant getArgumentData(int argument, Role role) override final;
   virtual void setArgument(int argument, QVariant value) override final;
   virtual void setArgument(int argument, QFile* file) override final;
   virtual void setArgument(int argument, EAbstractData* data) override final;
   quint32 getCapabilities() const override final { return Capabilities::Serial; }
   virtual bool initialize() override final;
   virtual void runSerial() override final;
   virtual int getBlockSize() override final { return 0; }
   virtual bool runBlock(int /*block*/) override final { return false; }
   virtual void finish() override final {}
private:
   using Transform = ExpressionMatrix::Transform;
   static const char* NONE;
   static const char* NLOG;
   static const char* LOG2;
   static const char* LOG10;
   QFile* _input {nullptr};
   ExpressionMatrix* _output {nullptr};
   QString _noSampleToken;
   qint32 _sampleSize {0};
   Transform _transform {Transform::None};
};



#endif
