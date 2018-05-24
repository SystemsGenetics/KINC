#ifndef EXPORTEXPRESSIONMATRIX_INPUT_H
#define EXPORTEXPRESSIONMATRIX_INPUT_H
#include "exportexpressionmatrix.h"



class ExportExpressionMatrix::Input : public EAbstractAnalytic::Input
{
   Q_OBJECT
public:
   enum Argument
   {
      InputData = 0
      ,OutputFile
      ,NoSampleToken
      ,Total
   };
   explicit Input(ExportExpressionMatrix* parent);
   virtual int size() const override final;
   virtual EAbstractAnalytic::Input::Type type(int index) const override final;
   virtual QVariant data(int index, Role role) const override final;
   virtual void set(int index, const QVariant& value) override final;
   virtual void set(int index, EAbstractData* data) override final;
   virtual void set(int index, QFile* file) override final;
private:
   ExportExpressionMatrix* _base;
};



#endif
