#ifndef EXPORTPARAMETERMATRIX_INPUT_H
#define EXPORTPARAMETERMATRIX_INPUT_H
#include "exportparametermatrix.h"



/*!
 * This class implements the abstract input of the export correlation matrix analytic.
 */
class ExportParameterMatrix::Input : public EAbstractAnalyticInput
{
   Q_OBJECT
public:
   /*!
    * Defines all input arguments for this analytic.
    */
   enum Argument
   {
      ExpressionData = 0
      ,ClusterData
      ,ParameterData
      ,Total
   };
   explicit Input(ExportParameterMatrix* parent);
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
   ExportParameterMatrix* _base;
};



#endif
