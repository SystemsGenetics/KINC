#ifndef EXPORTCORRELATIONMATRIX_INPUT_H
#define EXPORTCORRELATIONMATRIX_INPUT_H
#include "exportcorrelationmatrix.h"



/*!
 * This class implements the abstract input of the export correlation matrix analytic.
 */
class ExportCorrelationMatrix::Input : public EAbstractAnalytic::Input
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
      ,CorrelationData
      ,OutputFile
      ,Total
   };
   explicit Input(ExportCorrelationMatrix* parent);
   virtual int size() const override final;
   virtual EAbstractAnalytic::Input::Type type(int index) const override final;
   virtual QVariant data(int index, Role role) const override final;
   virtual void set(int index, const QVariant& value) override final;
   virtual void set(int index, EAbstractData* data) override final;
   virtual void set(int index, QFile* file) override final;
private:
   /*!
    * Pointer to the base analytic for this object.
    */
   ExportCorrelationMatrix* _base;
};



#endif
