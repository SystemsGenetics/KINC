#include "conditionaltest_input.h"
#include "datafactory.h"



/*!
 * Create a new input object.
 *
 * @param parent The parent analytic for this input object
 */
ConditionalTest::Input::Input(ConditionalTest* parent)
    : EAbstractAnalyticInput(parent), _base(parent)
{
    EDEBUG_FUNC(this,parent);
}



/*!
 * Query for the number of total inputs for the parent analytic.
 *
 * @return The integer representation for the total number of inputs.
 */
int ConditionalTest::Input::size() const
{
    EDEBUG_FUNC(this);

    return Total;
}



/*!
 * Query for the data types of the inputs.
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
    case EMXINPUT: return DataIn;
    case CCMINPUT: return DataIn;
    case CMXINPUT: return DataIn;
    case AMXINPUT: return FileIn;
    case Delimiter: return String;
    case MISSING: return String;
    case CSMOUT: return DataOut;
    case OVERRIDES: return String;
    case TEST: return String;
    default: return Boolean;
    }
}



/*!
 * Query infomation about the inputs.
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
    case EMXINPUT:
        switch(role)
        {
        case CommandLineName: return QString("emx");
        case Title: return tr("Input Expression Matrix:");
        case WhatsThis: return tr("A data file created by KINC containing the gene expression matrix created by the Import Expression Matrix analytic.");
        case DataType: return DataFactory::ExpressionMatrixType;
        default: return QVariant();
        }
    case CCMINPUT:
        switch(role)
        {
        case CommandLineName: return QString("ccm");
        case Title: return tr("Input Cluster Matrix:");
        case WhatsThis: return tr("A data file created by KINC containing the cluster sample masks created by the similarity analytic.");
        case DataType: return DataFactory::CCMatrixType;
        default: return QVariant();
        }
    case CMXINPUT:
        switch(role)
        {
        case CommandLineName: return QString("cmx");
        case Title: return tr("Input Correlation Matrix:");
        case WhatsThis: return tr("A data file created by KINC containing the correlation matrix values created by the similarity analytic.");
        case DataType: return DataFactory::CorrelationMatrixType;
        default: return QVariant();
        }
    case AMXINPUT:
        switch(role)
        {
        case CommandLineName: return QString("amx");
        case Title: return tr("Input Annotation matrix:");
        case WhatsThis: return tr("Tab delimited file where the row names represent samples in the experiment and the column names are conditional features of the experiment and the values are observed or measured values in the experiment.");
        case FileFilters: return tr("Annotation Matrix (*.txt)");
        default: return QVariant();
        }
    case Delimiter:
        switch(role)
        {
        case CommandLineName: return QString("delim");
        case Title: return tr("Annotation Matrix Delimiter:");
        case WhatsThis: return tr("Delimiter used to seperate values in the Annotation matrix, usually a tab or comma");
        case Default: return tr("tab");
        default: return QVariant();
        }
    case MISSING:
        switch(role)
        {
        case CommandLineName: return QString("missing");
        case Title: return tr("Missing value:");
        case WhatsThis: return tr("The string that specifies the missing value in the annotation matrix (e.g. NA, 0, 0.0).");
        case Default: return tr("NA");
        default: return QVariant();
        }
    case CSMOUT:
        EDEBUG_FUNC(this, role);
        switch(role)
        {
        case CommandLineName: return QString("output");
        case Title: return tr("Output Condition-Specific Matrix:");
        case WhatsThis: return tr("Condition-Specific Martrix, contains a matrix of clusters and their corrosponding p-values.");
        case DataType: return DataFactory::CSMatrixType;
        default: return QVariant();
        }
    case OVERRIDES:
        switch(role)
        {
        case CommandLineName: return QString("feat-types");
        case Title: return tr("Feature Types:");
        case WhatsThis: return tr("By default, this program will automatically detect the type of feature as 'categorical', 'quantitative', or 'ordinal'.   You can override the default type by listing the column name from the annotation matrix, followed by a colon and then the desired type. You can list as many features by separating them with commas, with no spaces around commas. For example if a column is named \"Health_Status\" and is numeric with an ordinal enter:  Health_Status:ordinal");
        default: return QVariant();
        }
    case TEST:
        switch(role)
        {
        case CommandLineName: return QString("feat-tests");
        case Title: return tr("Features to Test:");
        case WhatsThis: return tr("A comma-separated list of features, with no spaces around commas, from column names of the annotation matrix that should be tested. For example, if the annotation matrix has columns 'Treatment' and 'Subspecies' you can enter: \"Treatment,Subspecies\" Note: column names are case-sensitive.");
        default: return QVariant();
        }
    default: return QVariant();
    }
}



/*!
 * Set values to their respective places once they are inputed.
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
    case Delimiter:
        _base->_delimiter = value.toString();
        break;
    case MISSING:
        _base->_missing = value.toString();
        break;
    case TEST:
        _base->_userTestsStr = value.toString();
        break;
    case OVERRIDES:
        _base->_userTestTypesStr = value.toString();
        break;
    }
}



/*!
 * Set values to their respective places once they are inputed.
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
    case AMXINPUT:
        _base->_amx = file;
        break;
    }
}



/*!
 * Set values to their respective places once they are inputed.
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
    case EMXINPUT:
        _base->_emx = data->cast<ExpressionMatrix>();
        break;
    case CCMINPUT:
        _base->_ccm = data->cast<CCMatrix>();
        break;
    case CMXINPUT:
        _base->_cmx = data->cast<CorrelationMatrix>();
        break;
    case CSMOUT:
        _base->_out = data->cast<CSMatrix>();
        break;
    }
}
