#ifndef IMPORTEXPRESSIONMATRIX_H
#define IMPORTEXPRESSIONMATRIX_H
#include <ace/core/core.h>

#include "expressionmatrix.h"



class ImportExpressionMatrix : public EAbstractAnalytic
{
   Q_OBJECT
public:
   class Input;
   virtual int size() const override final;
   virtual void process(const EAbstractAnalytic::Block* result) override final;
   virtual EAbstractAnalytic::Input* makeInput() override final;
   virtual void initialize();
private:
   QFile* _input {nullptr};
   ExpressionMatrix* _output {nullptr};
   QString _noSampleToken;
   qint32 _sampleSize {0};
};



#endif
