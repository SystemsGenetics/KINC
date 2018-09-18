#ifndef IMPORTEXPRESSIONMATRIX_INPUT_H
#define IMPORTEXPRESSIONMATRIX_INPUT_H
#include "importexpressionmatrix.h"



class ImportExpressionMatrix::Input : public EAbstractAnalytic::Input
{
   Q_OBJECT
public:
   enum Argument
   {
      InputFile = 0
      ,OutputData
      ,NoSampleToken
      ,SampleSize
      ,Total
   };
   explicit Input(ImportExpressionMatrix* parent);
   virtual int size() const override final;
   virtual EAbstractAnalytic::Input::Type type(int index) const override final;
   virtual QVariant data(int index, Role role) const override final;
   virtual void set(int index, const QVariant& value) override final;
   virtual void set(int index, QFile* file) override final;
   virtual void set(int index, EAbstractData* data) override final;
private:
   ImportExpressionMatrix* _base;
};



#endif
