#ifndef IMPORTEXPRESSIONMATRIX_INPUT_H
#define IMPORTEXPRESSIONMATRIX_INPUT_H
#include "importexpressionmatrix.h"



/*!
 * This class implements the abstract input of the import expression matrix analytic.
 */
class ImportExpressionMatrix::Input : public EAbstractAnalyticInput
{
   Q_OBJECT
public:
   /*!
    * Defines all input arguments for this analytic.
    */
   enum Argument
   {
      InputFile = 0
      ,OutputData
      ,NANToken
      ,SampleSize
      ,Total
   };
   explicit Input(ImportExpressionMatrix* parent);
   virtual int size() const override final;
   virtual EAbstractAnalyticInput::Type type(int index) const override final;
   virtual QVariant data(int index, Role role) const override final;
   virtual void set(int index, const QVariant& value) override final;
   virtual void set(int index, QFile* file) override final;
   virtual void set(int index, EAbstractData* data) override final;
private:
   /*!
    * Pointer to the base analytic for this object.
    */
   ImportExpressionMatrix* _base;
};



#endif
