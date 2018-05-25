#ifndef EXTRACT_INPUT_H
#define EXTRACT_INPUT_H
#include "extract.h"



class Extract::Input : public EAbstractAnalytic::Input
{
   Q_OBJECT
public:
   enum Argument
   {
      ExpressionData = 0
      ,ClusterData
      ,CorrelationData
      ,OutputFile
      ,GraphMLFile
      ,MinCorrelation
      ,MaxCorrelation
      ,Total
   };
   explicit Input(Extract* parent);
   virtual int size() const override final;
   virtual EAbstractAnalytic::Input::Type type(int index) const override final;
   virtual QVariant data(int index, Role role) const override final;
   virtual void set(int index, const QVariant& value) override final;
   virtual void set(int index, EAbstractData* data) override final;
   virtual void set(int index, QFile* file) override final;
private:
   Extract* _base;
};



#endif
