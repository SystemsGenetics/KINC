#include "conditionaltest_input.h"
#include "datafactory.h"


/*!
*  Implements an interface to create a new input object.
*
* @param parent The parent analytic for this input object
*/
ConditionalTest::Input::Input(ConditionalTest* parent) : EAbstractAnalyticInput(parent), _base(parent)
{
    EDEBUG_FUNC(this,parent);
}






/*!
*  Implements an interface to quiery for the number of total inputs for the parent analytic.
*
* @return The integer representation for the total number of inputs.
*/
int ConditionalTest::Input::size() const
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
EAbstractAnalyticInput::Type ConditionalTest::Input::type(int index) const
{
    EDEBUG_FUNC(this, index);
    switch(index)
    {
    case EMXINPUT           : return DataIn;
    case CCMINPUT           : return DataIn;
    case CMXINPUT           : return DataIn;
    case ANXINPUT           : return FileIn;
    case CSCMOUT            : return DataOut;
    case ProbabilitySuccess : return Double;
    case OVERRIDES          : return String;
    case TEST               : return String;
    default                 : return Boolean;
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
QVariant ConditionalTest::Input::data(int index, Role role) const
{
    EDEBUG_FUNC(this, index, role);
    switch(index)
    {
    case EMXINPUT           : return emxData(role);
    case CCMINPUT           : return ccmData(role);
    case CMXINPUT           : return cmxData(role);
    case ANXINPUT           : return anxData(role);
    case CSCMOUT            : return CSCMData(role);
    case ProbabilitySuccess : return ProbabilitySuccessData(role);
    case OVERRIDES          : return overridesData(role);
    case TEST               : return testData(role);
    default                 : return QVariant();
    }
}





/*!
*  Implements an interface to set values to their respective places once they are inputed.
*
* @param index The input index you want to know the data type for.
*
* @param value The data inputed from the user.
*/
void ConditionalTest::Input::set(int index, const QVariant& value)
{
    EDEBUG_FUNC(this, value);
    switch(index)
    {
    case ProbabilitySuccess:
        _base->_probabilitySuccess = value.toDouble();
        break;
    case TEST:
        _base->_Testing = value.toString();
        break;
    case OVERRIDES:
        _base->_testOverride = value.toString();
        break;
    }
}





/*!
*  Implements an interface to set values to their respective places once they are inputed.
*
* @param index The input index you want to know the data type for.
*
* @param file The file inputed from the user.
*/
void ConditionalTest::Input::set(int index, QFile* file)
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
void ConditionalTest::Input::set(int index, EAbstractData* data)
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
QVariant ConditionalTest::Input::emxData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
    {
    case CommandLineName: return QString("emx");
    case Title          : return tr("Input emx Data File:");
    case WhatsThis      : return tr("A data file created by KINC containing \
                                    the gene expression matrix created by \
                                    the Import Expression Matrix analytic.");
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
QVariant ConditionalTest::Input::ccmData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
    {
    case CommandLineName: return QString("ccm");
    case Title          : return tr("Input ccm Data File:");
    case WhatsThis      : return tr("A data file created by KINC containing \
                                    the cluster sample masks created by \
                                    the similarity analytic.");
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
QVariant ConditionalTest::Input::cmxData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
    {
    case CommandLineName: return QString("cmx");
    case Title          : return tr("Input cmx Data File:");
    case WhatsThis      : return tr("A data file created by KINC containing \
                                     the correlation matrix values created by \
                                     the similarity analytic.");
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
QVariant ConditionalTest::Input::anxData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
    {
    case CommandLineName: return QString("amx");
    case Title          : return tr("Annotation matrix:");
    case WhatsThis      : return tr("Tab delimited file where the row names \
                                     represent samples in the experiment and \
                                     the column names are conditional features \
                                     of the experiment and the values are \
                                     observed or measured values in the experiment.");
    case FileFilters    : return tr("Annotation Matrix (*.txt)");
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
QVariant ConditionalTest::Input::CSCMData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
    {
    case CommandLineName: return QString("out");
    case Title          : return tr("Output Condition-Specific Matrix:");
    case WhatsThis      : return tr("Condition-Specific Martrix, contains \
                                     a matrix of clusters and their corrosponding\
                                     p-values.");
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
QVariant ConditionalTest::Input::ProbabilitySuccessData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
    {
    case CommandLineName: return QString("psuc");
    case Title          : return tr("Probability Success:");
    case WhatsThis      : return tr("For categorical features a binomial test \
                                     is performed. This value should contain the \
                                     probability of success for samples in a GMM \
                                     cluster being of the category label.");
    case Default        : return 0.75;
    case Minimum        : return 0;
    case Maximum        : return 1;
    default: return QVariant();
    }
}





/*!
*  Implements an interface to grab info about the test.
*
* @param role The role you are interested in knowing about.
*
* @return The information requested.
*/
QVariant ConditionalTest::Input::testData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
   {
    case CommandLineName: return QString("feat-tests");
    case Title          : return tr("Features to Test:");
    case WhatsThis      : return tr("A comma-separated list of features, with \
                                     no spaces around commas, from column \
                                     names of the annotation matrix that \
                                     should be tested. For example, if the \
                                     annotation matrix has columns 'Treatment' \
                                     and 'Subspecies' you can enter: \
                                     \"Treatment,Subspecies\" \
                                     Note: column names are case-sensitive.");
    default             : return QVariant();
    }
}





/*!
*  Implements an interface to grab info about the override infromation.
*
* @param role The role you are interested in knowing about.
*
* @return The information requested.
*/
QVariant ConditionalTest::Input::overridesData(Role role) const
{
    EDEBUG_FUNC(this, role);
    switch(role)
    {
    case CommandLineName: return QString("feat-types");
    case Title          : return tr("Feature Types:");
    case WhatsThis      : return tr("By default, this program will automatically \
                                     detect the type of feature as 'categorical', \
                                     'decimal', or 'ordinal'.   You can override \
                                     the default type by listing the column name \
                                     from the annotation matrix, two colons and \
                                     the desired type. You can list as many as \
                                     features by separating them with commas, \
                                     with no spaces around commas. \
                                     For example if a column is named \"Health_Status\" and is \
                                     numeric with an ordinal enter:  Health_Status::ordinal");
    default             : return QVariant();
    }
}
