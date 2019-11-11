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
   case ClusterDataIn:
      switch (role)
      {
      case Role::CommandLineName: return QString("ccm-in");
      case Role::Title: return tr("Input Cluster Matrix:");
      case Role::WhatsThis: return tr("A data file created by KINC containing the cluster sample masks created by the similarity analytic.");
      case Role::DataType: return DataFactory::CCMatrixType;
      default: return QVariant();
      }
   case CorrelationDataIn:
      switch (role)
      {
      case Role::CommandLineName: return QString("cmx-in");
      case Role::Title: return tr("Input Correlation Matrix:");
      case Role::WhatsThis: return tr("A data file created by KINC containing the correlation matrix values created by the similarity analytic.");
      case Role::DataType: return DataFactory::CorrelationMatrixType;
      default: return QVariant();
      }
   case ClusterDataOut:
      switch (role)
      {
      case Role::CommandLineName: return QString("ccm-out");
      case Role::Title: return tr("Output Cluster Matrix:");
      case Role::WhatsThis: return tr("A data file created by KINC containing the cluster sample masks created by the similarity analytic.");
      case Role::DataType: return DataFactory::CCMatrixType;
      default: return QVariant();
      }
   case CorrelationDataOut:
      switch (role)
      {
      case Role::CommandLineName: return QString("cmx-out");
      case Role::Title: return tr("Output Correlation Matrix:");
      case Role::WhatsThis: return tr("A data file created by KINC containing the correlation matrix values created by the similarity analytic.");
      case Role::DataType: return DataFactory::CorrelationMatrixType;
      default: return QVariant();
      }
   case PowerThresholdAlpha:
       switch (role)
       {
       case Role::CommandLineName: return QString("alpha");
       case Role::Title: return tr("Signficance Level");
       case Role::WhatsThis: return tr("The significance level (i.e. Type I error rate, alpha) for the power test.");
       case Role::Default: return 0.001;
       case Role::Minimum: return -std::numeric_limits<float>::infinity();
       case Role::Maximum: return +std::numeric_limits<float>::infinity();
       default: return QVariant();
       }
   case PowerThresholdPower:
       switch (role)
       {
       case Role::CommandLineName: return QString("power");
       case Role::Title: return tr("Power of test");
       case Role::WhatsThis: return tr("The power value (i.e. 1 minus Type II error rate, 1 minus beta) for the power test. For example, if the desired Type II error rate is 0.2, then this value should be 0.8.");
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

   if ( index == ClusterDataIn )
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
 * Set a file argument with the given index to the given qt file pointer. This
 * implementation does nothing because this analytic has no file arguments.
 *
 * @param index
 * @param file
 */
void CorrPowerFilter::Input::set(int, QFile*)
{
   EDEBUG_FUNC(this);
}
