#ifndef ConditionalTest_INPUT_H
#define ConditionalTest_INPUT_H
#include <ace/core/core.h>
#include "conditionaltest.h"



class ConditionalTest::Input : public EAbstractAnalyticInput
{
    Q_OBJECT
public:
    /*!
     * Defines all input arguments
     */
    enum Argument
    {
        EMXINPUT = 0
        ,CCMINPUT
        ,CMXINPUT
        ,AMXINPUT
        ,Delimiter
        ,MISSING
        ,CSMOUT
        ,TEST
        ,OVERRIDES
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
private:
    /*!
     * Pointer to the Inputs base parent object
     */
    ConditionalTest* _base;
};



#endif
