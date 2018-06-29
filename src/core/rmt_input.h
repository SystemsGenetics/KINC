#ifndef RMT_INPUT_H
#define RMT_INPUT_H
#include "rmt.h"



class RMT::Input : public EAbstractAnalytic::Input
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
      ,MinUnfoldingPace
      ,MaxUnfoldingPace
      ,HistogramBinSize
      ,Total
   };
   explicit Input(RMT* parent);
   virtual int size() const override final;
   virtual EAbstractAnalytic::Input::Type type(int index) const override final;
   virtual QVariant data(int index, Role role) const override final;
   virtual void set(int index, const QVariant& value) override final;
   virtual void set(int index, QFile* file) override final;
   virtual void set(int index, EAbstractData* data) override final;
private:
   RMT* _base;
};



#endif
