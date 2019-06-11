#include "corrpower_input.h"
#include "datafactory.h"



/*!
 * Construct a new input object with the given analytic as its parent.
 *
 * @param parent
 */
CorrPowerFilter::Input::Input(CorrPowerFilter* parent):
   EAbstractAnalyticInput(parent),
   _base(parent)
{
   EDEBUG_FUNC(this,parent);
}






/*!
 * Return the total number of arguments this analytic type contains.
 */
int CorrPowerFilter::Input::size() const
{
   EDEBUG_FUNC(this);

   return Total;
}






/*!
 * Return the argument type for a given index.
 *
 * @param index
 */
EAbstractAnalyticInput::Type CorrPowerFilter::Input::type(int index) const
{
   EDEBUG_FUNC(this,index);

   switch (index)
   {
   case ExpressionData: return Type::DataIn;
   case ClusterDataIn: return Type::DataIn;
   case CorrelationDataIn: return Type::DataIn;
   case ClusterDataOut: return Type::DataOut;
   case CorrelationDataOut: return Type::DataOut;
   case PowerThresholdAlpha: return Type::Double;
   case PowerThresholdPower: return Type::Double;
   default: return Type::Boolean;
   }
}






/*!
 * Return data for a given role on an argument with the given index.
 *
 * @param index
 * @param role
 */
QVariant CorrPowerFilter::Input::data(int index, Role role) const
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
   case ClusterDataIn:
      switch (role)
      {
      case Role::CommandLineName: return QString("ccm");
      case Role::Title: return tr("Input Cluster Matrix:");
      case Role::WhatsThis: return tr("Input cluster matrix containing cluster composition data.");
      case Role::DataType: return DataFactory::CCMatrixType;
      default: return QVariant();
      }
   case CorrelationDataIn:
      switch (role)
      {
      case Role::CommandLineName: return QString("cmx");
      case Role::Title: return tr("Input Correlation Matrix:");
      case Role::WhatsThis: return tr("Input correlation matrix containing correlation data.");
      case Role::DataType: return DataFactory::CorrelationMatrixType;
      default: return QVariant();
      }
   case ClusterDataOut:
      switch (role)
      {
      case Role::CommandLineName: return QString("ccmout");
      case Role::Title: return tr("Output Cluster Matrix:");
      case Role::WhatsThis: return tr("Output cluster matrix containing cluster composition data.");
      case Role::DataType: return DataFactory::CCMatrixType;
      default: return QVariant();
      }
   case CorrelationDataOut:
      switch (role)
      {
      case Role::CommandLineName: return QString("cmxout");
      case Role::Title: return tr("Output Correlation Matrix:");
      case Role::WhatsThis: return tr("Output correlation matrix containing correlation data.");
      case Role::DataType: return DataFactory::CorrelationMatrixType;
      default: return QVariant();
      }
   case PowerThresholdAlpha:
       switch (role)
       {
       case Role::CommandLineName: return QString("pwr.alpha");
       case Role::Title: return tr("Signficance Level (Type I error rate, alpha)");
       case Role::WhatsThis: return tr("If 'pwr.th' is TRUE then this is the Type I, alpha, significance level for the power analysis.");
       case Role::Default: return 0.001;
       case Role::Minimum: return -std::numeric_limits<float>::infinity();
       case Role::Maximum: return +std::numeric_limits<float>::infinity();
       default: return QVariant();
       }
   case PowerThresholdPower:
       switch (role)
       {
       case Role::CommandLineName: return QString("pwr.power");
       case Role::Title: return tr("The Power (1 - Type II error rate, 1-beta)");
       case Role::WhatsThis: return tr("If pwr.th is TRUE then this is the the power value (i.e. 1-Beta) for the power test.  If for example, the desired Type II error rate is 0.2, then this should be 0.8.");
       case Role::Default: return 0.8;
       case Role::Minimum: return -std::numeric_limits<float>::infinity();
       case Role::Maximum: return +std::numeric_limits<float>::infinity();
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
void CorrPowerFilter::Input::set(int index, const QVariant& value)
{
   EDEBUG_FUNC(this);
   switch (index)
   {
   case PowerThresholdAlpha:
      _base->_powerThresholdAlpha = value.toDouble();
      break;
   case PowerThresholdPower:
       _base->_powerThresholdPower = value.toDouble();
       break;
   }
}






/*!
 * Set a data argument with the given index to the given data object pointer.
 *
 * @param index
 * @param data
 */
void CorrPowerFilter::Input::set(int index, EAbstractData* data)
{
   EDEBUG_FUNC(this,index,data);

   if ( index == ExpressionData )
   {
      _base->_emx = data->cast<ExpressionMatrix>();
   }
   else if ( index == ClusterDataIn )
   {
      _base->_ccm = data->cast<CCMatrix>();
   }
   else if ( index == CorrelationDataIn )
   {
      _base->_cmx = data->cast<CorrelationMatrix>();
   }
   else if ( index ==  ClusterDataOut )
   {
      _base->_ccmOut = data->cast<CCMatrix>();
   }
   else if ( index == CorrelationDataOut )
   {
      _base->_cmxOut = data->cast<CorrelationMatrix>();
   }
}






/*!
 * Set a file argument with the given index to the given qt file pointer.
 *
 * @param index
 * @param file
 */
void CorrPowerFilter::Input::set(int index, QFile* file)
{
   EDEBUG_FUNC(this,index,file);

}
