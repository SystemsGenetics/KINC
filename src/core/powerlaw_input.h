#ifndef POWERLAW_INPUT_H
#define POWERLAW_INPUT_H
#include "powerlaw.h"



class PowerLaw::Input : public EAbstractAnalytic::Input
{
   Q_OBJECT
public:
   enum Argument
   {
      InputData = 0
      ,LogFile
      ,ThresholdStart
      ,ThresholdStep
      ,ThresholdStop
      ,Total
   };
   explicit Input(PowerLaw* parent);
   virtual int size() const override final;
   virtual EAbstractAnalytic::Input::Type type(int index) const override final;
   virtual QVariant data(int index, Role role) const override final;
   virtual void set(int index, const QVariant& value) override final;
   virtual void set(int index, QFile* file) override final;
   virtual void set(int index, EAbstractData* data) override final;
private:
   PowerLaw* _base;
};



#endif
