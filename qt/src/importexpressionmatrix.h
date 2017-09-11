#ifndef IMPORTEXPRESSIONMATRIX_H
#define IMPORTEXPRESSIONMATRIX_H
#include <ace/core/AceCore.h>



class ExpressionMatrix;



class ImportExpressionMatrix : public EAbstractAnalytic
{
public:
   enum Arguments
   {
      InputFile = 0
      ,OutputData
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
   virtual bool initialize() override final;
   virtual int getBlockSize() override final { return 1; }
   virtual bool runBlock(int block) override final;
   virtual void finish() override final {}
private:
   enum class Transform
   {
      None
      ,NLog
      ,Log2
      ,Log10
   };
   static const char* NONE;
   static const char* NLOG;
   static const char* LOG2;
   static const char* LOG10;
   QFile* _input {nullptr};
   ExpressionMatrix* _output {nullptr};
   QString _noSampleToken;//TODO: ADD ARGUMENT!
   qint32 _sampleSize {0};
   Transform _transform {Transform::None};
};



#endif
