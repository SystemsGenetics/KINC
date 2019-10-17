#include "extract_input.h"
#include "datafactory.h"



/*!
 * String list of output formats for this analytic that correspond exactly
 * to its enumeration. Used for handling the output format argument for this
 * input object.
 */
const QStringList Extract::Input::FORMAT_NAMES
{
   "text"
   ,"minimal"
   ,"graphml"
};






/*!
 * Construct a new input object with the given analytic as its parent.
 *
 * @param parent
 */
Extract::Input::Input(Extract* parent):
   EAbstractAnalyticInput(parent),
   _base(parent)
{
   EDEBUG_FUNC(this,parent);
}






/*!
 * Return the total number of arguments this analytic type contains.
 */
int Extract::Input::size() const
{
   EDEBUG_FUNC(this);

   return Total;
}






/*!
 * Return the argument type for a given index.
 *
 * @param index
 */
EAbstractAnalyticInput::Type Extract::Input::type(int index) const
{
   EDEBUG_FUNC(this,index);

   switch (index)
   {
   case ExpressionData: return Type::DataIn;
   case ClusterData: return Type::DataIn;
   case CorrelationData: return Type::DataIn;
   case ConditionSpecificClusterData : return Type::DataIn;
   case OutputFormatArg: return Type::Selection;
   case OutputFile: return Type::FileOut;
   case MinCorrelation: return Type::Double;
   case MaxCorrelation: return Type::Double;
   case CSMPValueFilter: return Type::String;
   default: return Type::Boolean;
   }
}






/*!
 * Return data for a given role on an argument with the given index.
 *
 * @param index
 * @param role
 */
QVariant Extract::Input::data(int index, Role role) const
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
   case ConditionSpecificClusterData :
       switch (role)
       {
       case Role::CommandLineName: return QString("csm");
       case Role::Title: return tr("Optional Condition Specific Cluster Matrix:");
      case Role::WhatsThis: return tr("Condition-Specific Martrix, contains a matrix of clusters and their corrosponding p-values.");
       case Role::DataType: return DataFactory::CSMatrixType;
       default: return QVariant();
       }
   case OutputFormatArg:
      switch (role)
      {
      case Role::CommandLineName: return QString("format");
      case Role::Title: return tr("Output Format:");
      case Role::WhatsThis: return tr("Format to use for the output file.");
      case Role::SelectionValues: return FORMAT_NAMES;
      case Role::Default: return "text";
      default: return QVariant();
      }
   case OutputFile:
      switch (role)
      {
      case Role::CommandLineName: return QString("output");
      case Role::Title: return tr("Output File:");
      case Role::WhatsThis: return tr("Output file that will contain network in the specified format.");
      default: return QVariant();
      }
   case MinCorrelation:
      switch (role)
      {
      case Role::CommandLineName: return QString("mincorr");
      case Role::Title: return tr("Minimum Correlation:");
      case Role::WhatsThis: return tr("Minimum (absolute) correlation threshold for gene pairs.");
      case Role::Default: return 0.85f;
      case Role::Minimum: return 0;
      case Role::Maximum: return 1;
      default: return QVariant();
      }
   case MaxCorrelation:
      switch (role)
      {
      case Role::CommandLineName: return QString("maxcorr");
      case Role::Title: return tr("Maximum Correlation:");
      case Role::WhatsThis: return tr("Maximum (absolute) correlation threshold for gene pairs.");
      case Role::Default: return 1;
      case Role::Minimum: return 0;
      case Role::Maximum: return 1;
      default: return QVariant();
      }
   case CSMPValueFilter:
      switch (role)
      {
      case Role::CommandLineName: return QString("filter-pvalue");
      case Role::Title: return tr("Optional P-Value Filter:");
      case Role::WhatsThis: return tr("An optional filter applied to the Condition-Specific Martrix provided above. The ouput network will not contain clusters with p-values above the given value for the given test labels. For example if you wanted to threshold at 1e-3 Subspecies Japonica, you would input \"Subspecies,Japonica,1e-3\", followed by \"::\" then any more p-value filters.");
      case Role::Default: return "";
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
void Extract::Input::set(int index, const QVariant& value)
{
   EDEBUG_FUNC(this,index,&value);

   switch (index)
   {
   case OutputFormatArg:
      _base->_outputFormat = static_cast<OutputFormat>(FORMAT_NAMES.indexOf(value.toString()));
      break;
   case MinCorrelation:
      _base->_minCorrelation = value.toFloat();
      break;
   case MaxCorrelation:
      _base->_maxCorrelation = value.toFloat();
       break;
    case CSMPValueFilter:
      _base->_csmPValueFilter = value.toString();
      break;
   }
}






/*!
 * Set a data argument with the given index to the given data object pointer.
 *
 * @param index
 * @param data
 */
void Extract::Input::set(int index, EAbstractData* data)
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
   else if ( index == ConditionSpecificClusterData )
   {
      _base->_csm = data->cast<CSMatrix>();
   }
}






/*!
 * Set a file argument with the given index to the given qt file pointer.
 *
 * @param index
 * @param file
 */
void Extract::Input::set(int index, QFile* file)
{
   EDEBUG_FUNC(this,index,file);

   if ( index == OutputFile )
   {
      _base->_output = file;
   }
}
