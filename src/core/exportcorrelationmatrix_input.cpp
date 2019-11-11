#include "exportcorrelationmatrix_input.h"
#include "datafactory.h"



/*!
 * Construct a new input object with the given analytic as its parent.
 *
 * @param parent
 */
ExportCorrelationMatrix::Input::Input(ExportCorrelationMatrix* parent):
   EAbstractAnalyticInput(parent),
   _base(parent)
{
   EDEBUG_FUNC(this,parent);
}



/*!
 * Return the total number of arguments this analytic type contains.
 */
int ExportCorrelationMatrix::Input::size() const
{
   EDEBUG_FUNC(this);

   return Total;
}



/*!
 * Return the argument type for a given index.
 *
 * @param index
 */
EAbstractAnalyticInput::Type ExportCorrelationMatrix::Input::type(int index) const
{
   EDEBUG_FUNC(this,index);

   switch (index)
   {
   case ExpressionData: return Type::DataIn;
   case ClusterData: return Type::DataIn;
   case CorrelationData: return Type::DataIn;
   case OutputFile: return Type::FileOut;
   default: return Type::Boolean;
   }
}



/*!
 * Return data for a given role on an argument with the given index.
 *
 * @param index
 * @param role
 */
QVariant ExportCorrelationMatrix::Input::data(int index, Role role) const
{
   EDEBUG_FUNC(this,index,role);

   switch (index)
   {
   case ExpressionData:
      switch (role)
      {
      case Role::CommandLineName: return QString("emx");
      case Role::Title: return tr("Input Expression Matrix:");
      case Role::WhatsThis: return tr("A data file created by KINC containing the gene expression matrix created by the Import Expression Matrix analytic.");
      case Role::DataType: return DataFactory::ExpressionMatrixType;
      default: return QVariant();
      }
   case ClusterData:
      switch (role)
      {
      case Role::CommandLineName: return QString("ccm");
      case Role::Title: return tr("Input Cluster Matrix:");
      case Role::WhatsThis: return tr("A data file created by KINC containing the cluster sample masks created by the similarity analytic.");
      case Role::DataType: return DataFactory::CCMatrixType;
      default: return QVariant();
      }
   case CorrelationData:
      switch (role)
      {
      case Role::CommandLineName: return QString("cmx");
      case Role::Title: return tr("Input Correlation Matrix:");
      case Role::WhatsThis: return tr("A data file created by KINC containing the correlation matrix values created by the similarity analytic.");
      case Role::DataType: return DataFactory::CorrelationMatrixType;
      default: return QVariant();
      }
   case OutputFile:
      switch (role)
      {
      case Role::CommandLineName: return QString("output");
      case Role::Title: return tr("Output File:");
      case Role::WhatsThis: return tr("Output text file that will contain pairwise correlation data.");
      case Role::FileFilters: return tr("Text file %1").arg("(*.txt)");
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
void ExportCorrelationMatrix::Input::set(int, const QVariant&)
{
   EDEBUG_FUNC(this);
}



/*!
 * Set a data argument with the given index to the given data object pointer.
 *
 * @param index
 * @param data
 */
void ExportCorrelationMatrix::Input::set(int index, EAbstractData* data)
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
   else if ( index == CorrelationData )
   {
      _base->_cmx = data->cast<CorrelationMatrix>();
   }
}



/*!
 * Set a file argument with the given index to the given qt file pointer.
 *
 * @param index
 * @param file
 */
void ExportCorrelationMatrix::Input::set(int index, QFile* file)
{
   EDEBUG_FUNC(this,index,file);

   if ( index == OutputFile )
   {
      _base->_output = file;
   }
}
