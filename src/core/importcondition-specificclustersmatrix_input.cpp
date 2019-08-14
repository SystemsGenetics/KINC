#include "importcondition-specificclustersmatrix_input.h"
#include "datafactory.h"


/*!
*  Implements an interface to create a new input object.
*
* @param parent The parent analytic for this input object
*/
importCSCM::Input::Input(importCSCM* parent) : EAbstractAnalyticInput(parent), _base(parent)
{
    EDEBUG_FUNC(this,parent);
}






/*!
*  Implements an interface to quiery for the number of total inputs for the parent analytic.
*
* @return The integer representation for the total number of inputs.
*/
int importCSCM::Input::size() const
{
    EDEBUG_FUNC(this);
    return Total;
}






/*!
*  Implements an interface to quiery for the data types of the inputs.
*
* @param index The input index you want to know the data type for.
*
* @return The type of the input.
*/
EAbstractAnalyticInput::Type importCSCM::Input::type(int index) const
{
    EDEBUG_FUNC(this, index);
    switch(index)
    {
    case EMXINPUT   : return DataIn;
    case CCMINPUT   : return DataIn;
    case CMXINPUT   : return DataIn;
    case ANXINPUT   : return FileIn;
    case CSCMOUT    : return DataOut;
    case ALPHA      : return Double;
    case OVERRIDES  : return String;
    case TEST       : return String;
    case CORRTHRESH : return Double;
    default         : return Boolean;
    }
}






/*!
*  Implements an interface to quiery infomation about the inputs.
*
* @param index The input index you want to know the data type for.
*
* @param role The role you are interested in knowing about.
*
* @return The information about the input in question.
*/
QVariant importCSCM::Input::data(int index, Role role) const
{
    EDEBUG_FUNC(this, index, role);
    switch(index)
    {
    case EMXINPUT   : return emxData(role);
    case CCMINPUT   : return ccmData(role);
    case CMXINPUT   : return cmxData(role);
    case ANXINPUT   : return anxData(role);
    case CSCMOUT    : return CSCMData(role);
    case ALPHA      : return alphaData(role);
    case OVERRIDES  : return overridesData(role);
    case TEST       : return testData(role);
    case CORRTHRESH : return corrThreshData(role);
    default         : return QVariant();
    }
}





/*!
*  Implements an interface to set values to their respective places once they are inputed.
*
* @param index The input index you want to know the data type for.
*
* @param value The data inputed from the user.
*/
void importCSCM::Input::set(int index, const QVariant& value)
{
    EDEBUG_FUNC(this, value);
    switch(index)
    {
    case ALPHA:
        _base->_alpha = value.toDouble();
        break;
    case TEST:
        _base->_Testing = value.toString();
        break;
    case OVERRIDES:
        _base->_testOverride = value.toString();
        break;
    case CORRTHRESH:
        _base->_corrthresh = value.toDouble();
    }
}





/*!
*  Implements an interface to set values to their respective places once they are inputed.
*
* @param index The input index you want to know the data type for.
*
* @param file The file inputed from the user.
*/
void importCSCM::Input::set(int index, QFile* file)
{
    EDEBUG_FUNC(this, index, file);
    switch(index)
    {
    case ANXINPUT:
        _base->_anx = file;
        break;
    }
}




/*!
*  Implements an interface to set values to their respective places once they are inputed.
*
* @param index The input index you want to know the data type for.
*
* @param value The data object inputed from the user.
*/
void importCSCM::Input::set(int index, EAbstractData* data)
{
    EDEBUG_FUNC(this, index, data);
    switch(index)
    {
    case EMXINPUT :
        _base->_emx = data->cast<ExpressionMatrix>();
        break;
    case CCMINPUT :
        _base->_ccm = data->cast<CCMatrix>();
        break;
    case CMXINPUT :
        _base->_cmx = data->cast<CorrelationMatrix>();
        break;
    case CSCMOUT:
        _base->_out = data->cast<CSCM>();
        break;
    }
}





/*!
*  Implements an interface to grab info about the emx data.
*
* @param role The role you are interested in knowing about.
*
* @return The information requested.
*/
QVariant importCSCM::Input::emxData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
    {
    case CommandLineName: return QString("emx");
    case Title          : return tr("Input emx Data File:");
    case WhatsThis      : return tr("emx table file");
    case DataType       : return DataFactory::ExpressionMatrixType;
    default             : return QVariant();
    }
}






/*!
*  Implements an interface to grab info about the ccm data.
*
* @param role The role you are interested in knowing about.
*
* @return The information requested.
*/
QVariant importCSCM::Input::ccmData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
    {
    case CommandLineName: return QString("ccm");
    case Title          : return tr("Input ccm Data File:");
    case WhatsThis      : return tr("ccm table file");
    case DataType       : return DataFactory::CCMatrixType;
    default             : return QVariant();
    }
}






/*!
*  Implements an interface to grab info about the cmx data.
*
* @param role The role you are interested in knowing about.
*
* @return The information requested.
*/
QVariant importCSCM::Input::cmxData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
    {
    case CommandLineName: return QString("cmx");
    case Title          : return tr("Input cmx Data File:");
    case WhatsThis      : return tr("cmx table file");
    case DataType       : return DataFactory::CorrelationMatrixType;
    default             : return QVariant();
    }
}






/*!
*  Implements an interface to grab info about the anx data.
*
* @param role The role you are interested in knowing about.
*
* @return The information requested.
*/
QVariant importCSCM::Input::anxData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
    {
    case CommandLineName: return QString("anx");
    case Title          : return tr("Input anx File:");
    case WhatsThis      : return tr("Raw anx file");
    case FileFilters    : return tr("Raw anx file (*.txt)");
    default             : return QVariant();
    }
}





/*!
*  Implements an interface to grab info about the CSCM data.
*
* @param role The role you are interested in knowing about.
*
* @return The information requested.
*/
QVariant importCSCM::Input::CSCMData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
    {
    case CommandLineName: return QString("out");
    case Title          : return tr("Output:");
    case WhatsThis      : return tr("CSCM, containis a matrix of clusters and their p-values");
    case DataType       : return DataFactory::CSCMType;
    default             : return QVariant();
    }
}






/*!
*  Implements an interface to grab info about the alpha value.
*
* @param role The role you are interested in knowing about.
*
* @return The information requested.
*/
QVariant importCSCM::Input::alphaData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
    {
    case CommandLineName: return QString("alpha");
    case Title          : return tr("Alpha:");
    case WhatsThis      : return tr("Threshold for keeping data");
    case Default        : return 0;
    case Minimum        : return 0;
    case Maximum        : return 1;
    default: return QVariant();
    }
}




/*!
*  Implements an interface to grab info about the override infromation.
*
* @param role The role you are interested in knowing about.
*
* @return The information requested.
*/
QVariant importCSCM::Input::overridesData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
    {
    case CommandLineName: return QString("override");
    case Title          : return tr("Test Overrides:");
    case WhatsThis      : return tr("Tests to override, taken as \"feature::test\",...");
    default             : return QVariant();
    }
}





/*!
*  Implements an interface to grab info about the test.
*
* @param role The role you are interested in knowing about.
*
* @return The information requested.
*/
QVariant importCSCM::Input::testData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
    {
    case CommandLineName: return QString("test");
    case Title          : return tr("Test:");
    case WhatsThis      : return tr("Features to test.");
    default             : return QVariant();
    }
}





/*!
*  Implements an interface to grab info about the corrolation threshold.
*
* @param role The role you are interested in knowing about.
*
* @return The information requested.
*/
QVariant importCSCM::Input::corrThreshData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
    {
    case CommandLineName: return QString("corrthresh");
    case Title          : return tr("Corrolation Threshold:");
    case WhatsThis      : return tr("Threshold for testing data");
    case Default        : return 0.85;
    case Minimum        : return 0;
    case Maximum        : return 1;
    default             : return QVariant();
    }
}
