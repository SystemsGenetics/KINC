#ifndef POWERLAW_INPUT_H
#define POWERLAW_INPUT_H
#include "powerlaw.h"



/*!
 * This class implements the abstract input of the PowerLaw analytic.
 */
class PowerLaw::Input : public EAbstractAnalyticInput
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
        ,ThresholdStart
        ,ThresholdStep
        ,ThresholdStop
        ,Total
    };
    explicit Input(PowerLaw* parent);
    virtual int size() const override final;
    virtual EAbstractAnalyticInput::Type type(int index) const override final;
    virtual QVariant data(int index, Role role) const override final;
    virtual void set(int index, const QVariant& value) override final;
    virtual void set(int index, QFile* file) override final;
    virtual void set(int index, EAbstractData* data) override final;
private:
    /*!
     * Pointer to the base analytic for this object.
     */
    PowerLaw* _base;
};



#endif
