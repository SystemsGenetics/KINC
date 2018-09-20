#include "rmt_input.h"
#include "correlationmatrix.h"
#include "datafactory.h"






/*!
 * Construct a new input object with the given analytic as its parent.
 *
 * @param parent
 */
RMT::Input::Input(RMT* parent):
   EAbstractAnalytic::Input(parent),
   _base(parent)
{}






/*!
 * Return the total number of arguments this analytic type contains.
 */
int RMT::Input::size() const
{
   return Total;
}






/*!
 * Return the argument type for a given index.
 *
 * @param index
 */
EAbstractAnalytic::Input::Type RMT::Input::type(int index) const
{
   switch (index)
   {
   case InputData: return Type::DataIn;
   case LogFile: return Type::FileOut;
   case ThresholdStart: return Type::Double;
   case ThresholdStep: return Type::Double;
   case ThresholdStop: return Type::Double;
   case MinUnfoldingPace: return Type::Integer;
   case MaxUnfoldingPace: return Type::Integer;
   case HistogramBinSize: return Type::Integer;
   default: return Type::Boolean;
   }
}






/*!
 * Return data for a given role on an argument with the given index.
 *
 * @param index
 * @param role
 */
QVariant RMT::Input::data(int index, Role role) const
{
   switch (index)
   {
   case InputData:
      switch (role)
      {
      case Role::CommandLineName: return QString("input");
      case Role::Title: return tr("Input:");
      case Role::WhatsThis: return tr("Correlation matrix for which an appropriate correlation threshold will be found.");
      case Role::DataType: return DataFactory::CorrelationMatrixType;
      default: return QVariant();
      }
   case LogFile:
      switch (role)
      {
      case Role::CommandLineName: return QString("log");
      case Role::Title: return tr("Log File:");
      case Role::WhatsThis: return tr("Output text file that logs all results.");
      case Role::FileFilters: return tr("Text file %1").arg("(*.txt)");
      default: return QVariant();
      }
   case ThresholdStart:
      switch (role)
      {
      case Role::CommandLineName: return QString("tstart");
      case Role::Title: return tr("Threshold Start:");
      case Role::WhatsThis: return tr("Starting threshold.");
      case Role::Default: return 0.99;
      case Role::Minimum: return 0;
      case Role::Maximum: return 1;
      default: return QVariant();
      }
   case ThresholdStep:
      switch (role)
      {
      case Role::CommandLineName: return QString("tstep");
      case Role::Title: return tr("Threshold Step:");
      case Role::WhatsThis: return tr("Threshold step size.");
      case Role::Default: return 0.001;
      case Role::Minimum: return 0;
      case Role::Maximum: return 1;
      default: return QVariant();
      }
   case ThresholdStop:
      switch (role)
      {
      case Role::CommandLineName: return QString("tstop");
      case Role::Title: return tr("Threshold Stop:");
      case Role::WhatsThis: return tr("Stopping threshold.");
      case Role::Default: return 0.5;
      case Role::Minimum: return 0;
      case Role::Maximum: return 1;
      default: return QVariant();
      }
   case MinUnfoldingPace:
      switch (role)
      {
      case Role::CommandLineName: return QString("minpace");
      case Role::Title: return tr("Minimum Unfolding Pace:");
      case Role::WhatsThis: return tr("The minimum pace with which to perform unfolding.");
      case Role::Default: return 10;
      case Role::Minimum: return 1;
      case Role::Maximum: return std::numeric_limits<int>::max();
      default: return QVariant();
      }
   case MaxUnfoldingPace:
      switch (role)
      {
      case Role::CommandLineName: return QString("maxpace");
      case Role::Title: return tr("Maximum Unfolding Pace:");
      case Role::WhatsThis: return tr("The maximum pace with which to perform unfolding.");
      case Role::Default: return 40;
      case Role::Minimum: return 1;
      case Role::Maximum: return std::numeric_limits<int>::max();
      default: return QVariant();
      }
   case HistogramBinSize:
      switch (role)
      {
      case Role::CommandLineName: return QString("bins");
      case Role::Title: return tr("Histogram Bin Size:");
      case Role::WhatsThis: return tr("The number of bins for the nearest-neighbor spacing histogram.");
      case Role::Default: return 60;
      case Role::Minimum: return 1;
      case Role::Maximum: return std::numeric_limits<int>::max();
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
void RMT::Input::set(int index, const QVariant& value)
{
   switch (index)
   {
   case ThresholdStart:
      _base->_thresholdStart = value.toDouble();
      break;
   case ThresholdStep:
      _base->_thresholdStep = value.toDouble();
      break;
   case ThresholdStop:
      _base->_thresholdStop = value.toDouble();
      break;
   case MinUnfoldingPace:
      _base->_minUnfoldingPace = value.toInt();
      break;
   case MaxUnfoldingPace:
      _base->_maxUnfoldingPace = value.toInt();
      break;
   case HistogramBinSize:
      _base->_histogramBinSize = value.toInt();
      break;
   }
}






/*!
 * Set a file argument with the given index to the given qt file pointer.
 *
 * @param index
 * @param file
 */
void RMT::Input::set(int index, QFile* file)
{
   if ( index == LogFile )
   {
      _base->_logfile = file;
   }
}






/*!
 * Set a data argument with the given index to the given data object pointer.
 *
 * @param index
 * @param data
 */
void RMT::Input::set(int index, EAbstractData* data)
{
   if ( index == InputData )
   {
      _base->_input = data->cast<CorrelationMatrix>();
   }
}
