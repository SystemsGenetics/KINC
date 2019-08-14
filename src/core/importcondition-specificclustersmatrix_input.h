#ifndef IMPORTCSCM_INPUT_H
#define IMPORTCSCM_INPUT_H
#include <ace/core/core.h>
#include "importcondition-specificclustersmatrix.h"

class importCSCM::Input : public EAbstractAnalyticInput
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
       *  defines CSCM modlie as an output
       */
      ,CSCMOUT
       /*!
       *  defines the threshold for data to keep
       */
      ,ALPHA
       /*!
       *  defines the corrolation threshold
       */
      ,CORRTHRESH
       /*!
       *  defines the overrides for the testtypes at a target location
       */
      ,OVERRIDES
       /*!
       *  defines the features not to test
       */
      ,TEST
       /*!
       *  defines Total number of arguments this class contains
       */
      ,Total
   };
public:
   explicit Input(importCSCM* parent);
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
   QVariant CSCMData(Role role) const;
   QVariant alphaData(Role role) const;
   QVariant overridesData(Role role) const;
   QVariant testData(Role role) const;
   QVariant corrThreshData(Role role) const;
private:
   /*!
   *  pointer to the Inputs base parent object
   */
   importCSCM* _base;
};

#endif // IMPORTCSCM_INPUT_H
