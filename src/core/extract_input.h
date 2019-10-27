#ifndef EXTRACT_INPUT_H
#define EXTRACT_INPUT_H
#include "extract.h"



/*!
 * This class implements the abstract input of the extract analytic.
 */
class Extract::Input : public EAbstractAnalyticInput
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
      ,ConditionSpecificClusterData
      ,OutputFormatArg
      ,OutputFile
      ,MinCorrelation
      ,MaxCorrelation
      ,CSMPValueFilter
      ,CSMRSquareFilter
      ,Total
   };
   explicit Input(Extract* parent);
   virtual int size() const override final;
   virtual EAbstractAnalyticInput::Type type(int index) const override final;
   virtual QVariant data(int index, Role role) const override final;
   virtual void set(int index, const QVariant& value) override final;
   virtual void set(int index, EAbstractData* data) override final;
   virtual void set(int index, QFile* file) override final;
private:
   static const QStringList FORMAT_NAMES;
   /*!
    * Pointer to the base analytic for this object.
    */
   Extract* _base;
};



#endif
