#ifndef IMPORTCORRELATIONMATRIX_INPUT_H
#define IMPORTCORRELATIONMATRIX_INPUT_H
#include "importcorrelationmatrix.h"



class ImportCorrelationMatrix::Input : public EAbstractAnalytic::Input
{
   Q_OBJECT
public:
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
   ImportCorrelationMatrix* _base;
};



#endif
