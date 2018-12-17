#ifndef IMPORTCORRELATIONMATRIX_INPUT_H
#define IMPORTCORRELATIONMATRIX_INPUT_H
#include "importcorrelationmatrix.h"



/*!
 * This class implements the abstract input of the import correlation matrix analytic.
 */
class ImportCorrelationMatrix::Input : public EAbstractAnalytic::Input
{
   Q_OBJECT
public:
   /*!
    * Defines all input arguments for this analytic.
    */
   enum Argument
   {
      InputFile = 0
      ,ClusterData
      ,CorrelationData
      ,GeneSize
      ,MaxClusterSize
      ,SampleSize
      ,CorrelationName
      ,Total
   };
   explicit Input(ImportCorrelationMatrix* parent);
   virtual int size() const override final;
   virtual EAbstractAnalytic::Input::Type type(int index) const override final;
   virtual QVariant data(int index, Role role) const override final;
   virtual void set(int index, const QVariant& value) override final;
   virtual void set(int index, QFile* file) override final;
   virtual void set(int index, EAbstractData* data) override final;
private:
   /*!
    * Pointer to the base analytic for this object.
    */
   ImportCorrelationMatrix* _base;
};



#endif
