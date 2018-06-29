#ifndef EXPORTCORRELATIONMATRIX_INPUT_H
#define EXPORTCORRELATIONMATRIX_INPUT_H
#include "exportcorrelationmatrix.h"



class ExportCorrelationMatrix::Input : public EAbstractAnalytic::Input
{
   Q_OBJECT
public:
   enum Argument
   {
      ClusterData = 0
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
   ExportCorrelationMatrix* _base;
};



#endif
