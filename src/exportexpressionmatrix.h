#ifndef EXPORTEXPRESSIONMATRIX_H
#define EXPORTEXPRESSIONMATRIX_H
#include <ace/core/core.h>

#include "expressionmatrix.h"



class ExportExpressionMatrix : public EAbstractAnalytic
{
   Q_OBJECT
public:
   class Input;
   virtual int size() const override final;
   virtual void process(const EAbstractAnalytic::Block* result) override final;
   virtual EAbstractAnalytic::Input* makeInput() override final;
   virtual void initialize();
private:
   ExpressionMatrix* _input {nullptr};
   QFile* _output {nullptr};
   QString _noSampleToken;
};



#endif
