#include "powerlaw_input.h"
#include "correlationmatrix.h"
#include "datafactory.h"



/*!
 * Construct a new input object with the given analytic as its parent.
 *
 * @param parent
 */
PowerLaw::Input::Input(PowerLaw* parent):
    EAbstractAnalyticInput(parent),
    _base(parent)
{
    EDEBUG_FUNC(this,parent);
}



/*!
 * Return the total number of arguments this analytic type contains.
 */
int PowerLaw::Input::size() const
{
    EDEBUG_FUNC(this);

    return Total;
}



/*!
 * Return the argument type for a given index.
 *
 * @param index
 */
EAbstractAnalyticInput::Type PowerLaw::Input::type(int index) const
{
    EDEBUG_FUNC(this,index);

    switch (index)
    {
    case InputData: return Type::DataIn;
    case LogFile: return Type::FileOut;
    case ThresholdStart: return Type::Double;
    case ThresholdStep: return Type::Double;
    case ThresholdStop: return Type::Double;
    default: return Type::Boolean;
    }
}



/*!
 * Return data for a given role on an argument with the given index.
 *
 * @param index
 * @param role
 */
QVariant PowerLaw::Input::data(int index, Role role) const
{
    EDEBUG_FUNC(this,index,role);

    switch (index)
    {
    case InputData:
        switch (role)
        {
        case Role::CommandLineName: return QString("input");
        case Role::Title: return tr("Input:");
        case Role::WhatsThis: return tr("Correlation matrix for which an appropriate correlation threshold will be found.");
        case Role::DataType: return DataFactory::CorrelationMatrixType;
        default: return QVariant();
        }
    case LogFile:
        switch (role)
        {
        case Role::CommandLineName: return QString("log");
        case Role::Title: return tr("Log File:");
        case Role::WhatsThis: return tr("Output text file that logs all results.");
        case Role::FileFilters: return tr("Text file %1").arg("(*.txt)");
        default: return QVariant();
        }
    case ThresholdStart:
        switch (role)
        {
        case Role::CommandLineName: return QString("tstart");
        case Role::Title: return tr("Threshold Start:");
        case Role::WhatsThis: return tr("Starting threshold.");
        case Role::Default: return 0.99;
        case Role::Minimum: return 0;
        case Role::Maximum: return 1;
        default: return QVariant();
        }
    case ThresholdStep:
        switch (role)
        {
        case Role::CommandLineName: return QString("tstep");
        case Role::Title: return tr("Threshold Step:");
        case Role::WhatsThis: return tr("Threshold step size.");
        case Role::Default: return 0.01;
        case Role::Minimum: return 0;
        case Role::Maximum: return 1;
        default: return QVariant();
        }
    case ThresholdStop:
        switch (role)
        {
        case Role::CommandLineName: return QString("tstop");
        case Role::Title: return tr("Threshold Stop:");
        case Role::WhatsThis: return tr("Stopping threshold.");
        case Role::Default: return 0.5;
        case Role::Minimum: return 0;
        case Role::Maximum: return 1;
        default: return QVariant();
        }
    default: return QVariant();
    }
}



/*!
 * Set an argument with the given index to the given value.
 *
 * @param index
 * @param value
 */
void PowerLaw::Input::set(int index, const QVariant& value)
{
    EDEBUG_FUNC(this,index,&value);

    switch (index)
    {
    case ThresholdStart:
        _base->_thresholdStart = value.toFloat();
        break;
    case ThresholdStep:
        _base->_thresholdStep = value.toFloat();
        break;
    case ThresholdStop:
        _base->_thresholdStop = value.toFloat();
        break;
    }
}



/*!
 * Set a file argument with the given index to the given qt file pointer.
 *
 * @param index
 * @param file
 */
void PowerLaw::Input::set(int index, QFile* file)
{
    EDEBUG_FUNC(this,index,file);

    if ( index == LogFile )
    {
        _base->_logfile = file;
    }
}



/*!
 * Set a data argument with the given index to the given data object pointer.
 *
 * @param index
 * @param data
 */
void PowerLaw::Input::set(int index, EAbstractData* data)
{
    EDEBUG_FUNC(this,index,data);

    if ( index == InputData )
    {
        _base->_input = data->cast<CorrelationMatrix>();
    }
}
