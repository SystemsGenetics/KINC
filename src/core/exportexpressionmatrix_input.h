#ifndef EXPORTEXPRESSIONMATRIX_INPUT_H
#define EXPORTEXPRESSIONMATRIX_INPUT_H
#include "exportexpressionmatrix.h"



/*!
 * This class implements the abstract input of the export expression matrix analytic.
 */
class ExportExpressionMatrix::Input : public EAbstractAnalyticInput
{
   Q_OBJECT
public:
   /*!
    * Defines all input arguments for this analytic.
    */
   enum Argument
   {
      InputData = 0
      ,OutputFile
      ,NANToken
      ,Precision
      ,Total
   };
   explicit Input(ExportExpressionMatrix* parent);
   virtual int size() const override final;
   virtual EAbstractAnalyticInput::Type type(int index) const override final;
   virtual QVariant data(int index, Role role) const override final;
   virtual void set(int index, const QVariant& value) override final;
   virtual void set(int index, EAbstractData* data) override final;
   virtual void set(int index, QFile* file) override final;
private:
   /*!
    * Pointer to the base analytic for this object.
    */
   ExportExpressionMatrix* _base;
};



#endif
