#include "importcorrelationmatrix_input.h"
#include "datafactory.h"
#include "pairwise_index.h"



/*!
 * Construct a new input object with the given analytic as its parent.
 *
 * @param parent
 */
ImportCorrelationMatrix::Input::Input(ImportCorrelationMatrix* parent):
   EAbstractAnalytic::Input(parent),
   _base(parent)
{
   EDEBUG_FUNC(this,parent);
}






/*!
 * Return the total number of arguments this analytic type contains.
 */
int ImportCorrelationMatrix::Input::size() const
{
   EDEBUG_FUNC(this);

   return Total;
}






/*!
 * Return the argument type for a given index.
 *
 * @param index
 */
EAbstractAnalytic::Input::Type ImportCorrelationMatrix::Input::type(int index) const
{
   EDEBUG_FUNC(this,index);

   switch (index)
   {
   case InputFile: return Type::FileIn;
   case ClusterData: return Type::DataOut;
   case CorrelationData: return Type::DataOut;
   case GeneSize: return Type::Integer;
   case MaxClusterSize: return Type::Integer;
   case SampleSize: return Type::Integer;
   case CorrelationName: return Type::String;
   default: return Type::Boolean;
   }
}






/*!
 * Return data for a given role on an argument with the given index.
 *
 * @param index
 * @param role
 */
QVariant ImportCorrelationMatrix::Input::data(int index, Role role) const
{
   EDEBUG_FUNC(this,index,role);

   switch (index)
   {
   case InputFile:
      switch (role)
      {
      case Role::CommandLineName: return QString("input");
      case Role::Title: return tr("Input File:");
      case Role::WhatsThis: return tr("Input text file containing pairwise correlation data.");
      case Role::FileFilters: return tr("Text file %1").arg("(*.txt)");
      default: return QVariant();
      }
   case ClusterData:
      switch (role)
      {
      case Role::CommandLineName: return QString("ccm");
      case Role::Title: return tr("Cluster Matrix:");
      case Role::WhatsThis: return tr("Output cluster matrix that will contain cluster composition data.");
      case Role::DataType: return DataFactory::CCMatrixType;
      default: return QVariant();
      }
   case CorrelationData:
      switch (role)
      {
      case Role::CommandLineName: return QString("cmx");
      case Role::Title: return tr("Correlation Matrix:");
      case Role::WhatsThis: return tr("Output correlation matrix that will contain correlation data.");
      case Role::DataType: return DataFactory::CorrelationMatrixType;
      default: return QVariant();
      }
   case GeneSize:
      switch (role)
      {
      case Role::CommandLineName: return QString("genes");
      case Role::Title: return tr("Gene Size:");
      case Role::WhatsThis: return tr("Number of genes.");
      case Role::Minimum: return 1;
      case Role::Maximum: return std::numeric_limits<int>::max();
      default: return QVariant();
      }
   case MaxClusterSize:
      switch (role)
      {
      case Role::CommandLineName: return QString("maxclusters");
      case Role::Title: return tr("Maximum Cluster Size:");
      case Role::WhatsThis: return tr("Maximum number of clusters per pair.");
      case Role::Minimum: return 1;
      case Role::Maximum: return Pairwise::Index::MAX_CLUSTER_SIZE;
      default: return QVariant();
      }
   case SampleSize:
      switch (role)
      {
      case Role::CommandLineName: return QString("samples");
      case Role::Title: return tr("Sample Size:");
      case Role::WhatsThis: return tr("Number of samples.");
      case Role::Minimum: return 1;
      case Role::Maximum: return std::numeric_limits<int>::max();
      default: return QVariant();
      }
   case CorrelationName:
      switch (role)
      {
      case Role::CommandLineName: return QString("corrname");
      case Role::Title: return tr("Correlation Name:");
      case Role::WhatsThis: return tr("Name of correlation method.");
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
void ImportCorrelationMatrix::Input::set(int index, const QVariant& value)
{
   EDEBUG_FUNC(this,index,value);

   switch (index)
   {
   case GeneSize:
      _base->_geneSize = value.toInt();
      break;
   case MaxClusterSize:
      _base->_maxClusterSize = value.toInt();
      break;
   case SampleSize:
      _base->_sampleSize = value.toInt();
      break;
   case CorrelationName:
      _base->_correlationName = value.toString();
      break;
   }
}






/*!
 * Set a file argument with the given index to the given qt file pointer.
 *
 * @param index
 * @param file
 */
void ImportCorrelationMatrix::Input::set(int index, QFile* file)
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
void ImportCorrelationMatrix::Input::set(int index, EAbstractData* data)
{
   EDEBUG_FUNC(this,index,data);

   if ( index == ClusterData )
   {
      _base->_ccm = data->cast<CCMatrix>();
   }
   else if ( index == CorrelationData )
   {
      _base->_cmx = data->cast<CorrelationMatrix>();
   }
}
