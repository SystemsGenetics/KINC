#include "exportexpressionmatrix_input.h"
#include "datafactory.h"






/*!
 * Construct a new input object with the given analytic as its parent.
 *
 * @param parent
 */
ExportExpressionMatrix::Input::Input(ExportExpressionMatrix* parent):
   EAbstractAnalytic::Input(parent),
   _base(parent)
{}






/*!
 * Return the total number of arguments this analytic type contains.
 */
int ExportExpressionMatrix::Input::size() const
{
   return Total;
}






/*!
 * Return the argument type for a given index.
 *
 * @param index
 */
EAbstractAnalytic::Input::Type ExportExpressionMatrix::Input::type(int index) const
{
   switch (index)
   {
   case InputData: return Type::DataIn;
   case OutputFile: return Type::FileOut;
   case NANToken: return Type::String;
   default: return Type::Boolean;
   }
}






/*!
 * Return data for a given role on an argument with the given index.
 *
 * @param index
 * @param role
 */
QVariant ExportExpressionMatrix::Input::data(int index, Role role) const
{
   switch (index)
   {
   case InputData:
      switch (role)
      {
      case Role::CommandLineName: return QString("input");
      case Role::Title: return tr("Input:");
      case Role::WhatsThis: return tr("Input expression matrix containing expression data.");
      case Role::DataType: return DataFactory::ExpressionMatrixType;
      default: return QVariant();
      }
   case OutputFile:
      switch (role)
      {
      case Role::CommandLineName: return QString("output");
      case Role::Title: return tr("Output:");
      case Role::WhatsThis: return tr("Output text file that will contain space/tab divided gene expression data.");
      case Role::FileFilters: return tr("Text file %1").arg("(*.txt)");
      default: return QVariant();
      }
   case NANToken:
      switch (role)
      {
      case Role::CommandLineName: return QString("nan");
      case Role::Title: return tr("NAN Token:");
      case Role::WhatsThis: return tr("Expected token for expressions that have no value.");
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
void ExportExpressionMatrix::Input::set(int index, const QVariant& value)
{
   switch (index)
   {
   case NANToken:
      _base->_nanToken = value.toString();
      break;
   }
}






/*!
 * Set a data argument with the given index to the given data object pointer.
 *
 * @param index
 * @param data
 */
void ExportExpressionMatrix::Input::set(int index, EAbstractData* data)
{
   if ( index == InputData )
   {
      _base->_input = data->cast<ExpressionMatrix>();
   }
}






/*!
 * Set a file argument with the given index to the given qt file pointer.
 *
 * @param index
 * @param file
 */
void ExportExpressionMatrix::Input::set(int index, QFile* file)
{
   if ( index == OutputFile )
   {
      _base->_output = file;
   }
}
