#ifndef RMT_INPUT_H
#define RMT_INPUT_H
#include "rmt.h"



/*!
 * This class implements the abstract input of the RMT analytic.
 */
class RMT::Input : public EAbstractAnalytic::Input
{
   Q_OBJECT
public:
   /*!
    * Defines all input arguments for this analytic.
    */
   enum Argument
   {
      InputData = 0
      ,LogFile
      ,ReductionType
      ,ThresholdStart
      ,ThresholdStep
      ,ThresholdStop
      ,SplineInterpolation
      ,MinSplinePace
      ,MaxSplinePace
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
   static const QStringList REDUCTION_NAMES;
   /*!
    * Pointer to the base analytic for this object.
    */
   RMT* _base;
};



#endif
