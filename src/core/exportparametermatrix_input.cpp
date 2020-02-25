#include "exportparametermatrix_input.h"
#include "datafactory.h"



/*!
 * Construct a new input object with the given analytic as its parent.
 *
 * @param parent
 */
ExportParameterMatrix::Input::Input(ExportParameterMatrix* parent):
    EAbstractAnalyticInput(parent),
    _base(parent)
{
    EDEBUG_FUNC(this,parent);
}






/*!
 * Return the total number of arguments this analytic type contains.
 */
int ExportParameterMatrix::Input::size() const
{
    EDEBUG_FUNC(this);

    return Total;
}






/*!
 * Return the argument type for a given index.
 *
 * @param index
 */
EAbstractAnalyticInput::Type ExportParameterMatrix::Input::type(int index) const
{
    EDEBUG_FUNC(this,index);

    switch (index)
    {
    case ExpressionData: return Type::DataIn;
    case ClusterData: return Type::DataIn;
    case ParameterData: return Type::DataOut;
    default: return Type::Boolean;
    }
}






/*!
 * Return data for a given role on an argument with the given index.
 *
 * @param index
 * @param role
 */
QVariant ExportParameterMatrix::Input::data(int index, Role role) const
{
    EDEBUG_FUNC(this,index,role);

    switch (index)
    {
    case ExpressionData:
        switch (role)
        {
        case Role::CommandLineName: return QString("emx");
        case Role::Title: return tr("Expression Matrix:");
        case Role::WhatsThis: return tr("Input expression matrix containing gene expression data.");
        case Role::DataType: return DataFactory::ExpressionMatrixType;
        default: return QVariant();
        }
    case ClusterData:
        switch (role)
        {
        case Role::CommandLineName: return QString("ccm");
        case Role::Title: return tr("Cluster Matrix:");
        case Role::WhatsThis: return tr("Input cluster matrix containing cluster composition data.");
        case Role::DataType: return DataFactory::CCMatrixType;
        default: return QVariant();
        }
    case ParameterData:
        switch (role)
        {
        case Role::CommandLineName: return QString("cpm");
        case Role::Title: return tr("Parameter Matrix:");
        case Role::WhatsThis: return tr("Output parameter matrix containing cluster parameter data.");
        case Role::DataType: return DataFactory::CPMatrixType;
        default: return QVariant();
        }
    default: return QVariant();
    }
}






/*!
 * Set an argument with the given index to the given value. This analytic has
 * no basic arguments so this function does nothing.
 *
 * @param index
 * @param value
 */
void ExportParameterMatrix::Input::set(int, const QVariant&)
{
    EDEBUG_FUNC(this);
}






/*!
 * Set a data argument with the given index to the given data object pointer.
 *
 * @param index
 * @param data
 */
void ExportParameterMatrix::Input::set(int index, EAbstractData* data)
{
    EDEBUG_FUNC(this,index,data);

    if ( index == ExpressionData )
    {
        _base->_emx = data->cast<ExpressionMatrix>();
    }
    else if ( index == ClusterData )
    {
        _base->_ccm = data->cast<CCMatrix>();
    }
    else if ( index == ParameterData )
    {
        _base->_cpm = data->cast<CPMatrix>();
    }
}






/*!
 * Set a file argument with the given index to the given qt file pointer.
 *
 * @param index
 * @param file
 */
void ExportParameterMatrix::Input::set(int, QFile*)
{
    EDEBUG_FUNC(this);
}
