#ifndef ConditionalTest_INPUT_H
#define ConditionalTest_INPUT_H
#include <ace/core/core.h>
#include "conditionaltest.h"

class ConditionalTest::Input : public EAbstractAnalyticInput
{
   Q_OBJECT
public:
    /*!
    *  Defines all arguments
    */
   enum Argument
   {
       /*!
       *  defines Input file argument for expresion matrix
       */
      EMXINPUT = 0
       /*!
       *  defines Input file argument for cluster composition matrix
       */
      ,CCMINPUT
       /*!
       *  defines Input file argument for correlation matrix
       */
      ,CMXINPUT
       /*!
       *  defines Input file argument for annotation matrix
       */
      ,ANXINPUT
       /*!
       *  defines CSM modlie as an output
       */
      ,CSMOUT
       /*!
       *  defines the threshold for data to keep
       */
      ,ProbabilitySuccess
       /*!
       *  defines the features not to test
       */
      ,TEST
       /*!
       *  defines the overrides for the testtypes at a target location
       */
      ,OVERRIDES
       /*!
       *  defines Total number of arguments this class contains
       */
      ,Total
   };
public:
   explicit Input(ConditionalTest* parent);
   virtual int size() const override final;
   virtual EAbstractAnalyticInput::Type type(int index) const override final;
   virtual QVariant data(int index, Role role) const override final;
   virtual void set(int index, const QVariant& value) override final;
   virtual void set(int index, QFile* file) override final;
   virtual void set(int index, EAbstractData* data) override final;
   QVariant emxData(Role role) const;
   QVariant ccmData(Role role) const;
   QVariant cmxData(Role role) const;
   QVariant anxData(Role role) const;
   QVariant CSMData(Role role) const;
   QVariant ProbabilitySuccessData(Role role) const;
   QVariant testData(Role role) const;
   QVariant overridesData(Role role) const;

private:
   /*!
   *  pointer to the Inputs base parent object
   */
   ConditionalTest* _base;
};

#endif // ConditionalTest_INPUT_H
