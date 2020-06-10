#include "importexpressionmatrix_input.h"
#include "datafactory.h"



/*!
 * Construct a new input object with the given analytic as its parent.
 *
 * @param parent
 */
ImportExpressionMatrix::Input::Input(ImportExpressionMatrix* parent):
    EAbstractAnalyticInput(parent),
    _base(parent)
{
    EDEBUG_FUNC(this,parent);
}



/*!
 * Return the total number of arguments this analytic type contains.
 */
int ImportExpressionMatrix::Input::size() const
{
    EDEBUG_FUNC(this);

    return Total;
}



/*!
 * Return the argument type for a given index.
 *
 * @param index
 */
EAbstractAnalyticInput::Type ImportExpressionMatrix::Input::type(int index) const
{
    EDEBUG_FUNC(this,index);

    switch (index)
    {
    case InputFile: return Type::FileIn;
    case OutputData: return Type::DataOut;
    case NANToken: return Type::String;
    case SampleSize: return Type::Integer;
    case ContainsRowID: return Type::Boolean;
    default: return Type::Boolean;
    }
}



/*!
 * Return data for a given role on an argument with the given index.
 *
 * @param index
 * @param role
 */
QVariant ImportExpressionMatrix::Input::data(int index, Role role) const
{
    EDEBUG_FUNC(this,index,role);

    switch (index)
    {
    case InputFile:
        switch (role)
        {
        case Role::CommandLineName: return QString("input");
        case Role::Title: return tr("Input:");
        case Role::WhatsThis: return tr("Input text file containing space/tab delimited gene expression data.");
        case Role::FileFilters: return tr("Text file %1").arg("(*.txt)");
        default: return QVariant();
        }
    case OutputData:
        switch (role)
        {
        case Role::CommandLineName: return QString("output");
        case Role::Title: return tr("Output Expression Matrix:");
        case Role::WhatsThis: return tr("A data file created by KINC containing the gene expression matrix created by the Import Expression Matrix analytic.");
        case Role::DataType: return DataFactory::ExpressionMatrixType;
        default: return QVariant();
        }
    case NANToken:
        switch (role)
        {
        case Role::CommandLineName: return QString("nan");
        case Role::Title: return tr("NAN Token:");
        case Role::WhatsThis: return tr("The string that specifies the missing value in the annotation matrix (e.g. NA, 0, 0.0).");
        case Role::Default: return "NA";
        default: return QVariant();
        }
    case SampleSize:
        switch (role)
        {
        case Role::CommandLineName: return QString("samples");
        case Role::Title: return tr("Sample Size:");
        case Role::WhatsThis: return tr("Number of samples. 0 indicates the text file contains a header of sample names to be read to determine size.");
        case Role::Default: return 0;
        case Role::Minimum: return 0;
        case Role::Maximum: return std::numeric_limits<int>::max();
        default: return QVariant();
        }
    case ContainsRowID:
        switch (role)
        {
        case Role::CommandLineName: return QString("contains-row-id");
        case Role::Title: return tr("Contains Row ID:");
        case Role::WhatsThis: return tr("Indicates whether the input expression matrix has a \"RowID\" entry that should be ignored.");
        case Role::Default: return false;
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
void ImportExpressionMatrix::Input::set(int index, const QVariant& value)
{
    EDEBUG_FUNC(this,index,&value);

    switch (index)
    {
    case NANToken:
        _base->_nanToken = value.toString();
        break;
    case SampleSize:
        _base->_sampleSize = value.toInt();
        break;
    case ContainsRowID:
        _base->_containsRowID = value.toBool();
        break;
    }
}



/*!
 * Set a file argument with the given index to the given qt file pointer.
 *
 * @param index
 * @param file
 */
void ImportExpressionMatrix::Input::set(int index, QFile* file)
{
    EDEBUG_FUNC(this,index,file);

    if ( index == InputFile )
    {
        _base->_input = file;
    }
}



/*!
 * Set a data argument with the given index to the given data object pointer.
 *
 * @param index
 * @param data
 */
void ImportExpressionMatrix::Input::set(int index, EAbstractData* data)
{
    EDEBUG_FUNC(this,index,data);

    if ( index == OutputData )
    {
        _base->_output = data->cast<ExpressionMatrix>();
    }
}
