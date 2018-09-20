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
{
   EDEBUG_FUNC(this,parent);
}






/*!
 * Return the total number of arguments this analytic type contains.
 */
int RMT::Input::size() const
{
   EDEBUG_FUNC(this);

   return Total;
}






/*!
 * Return the argument type for a given index.
 *
 * @param index
 */
EAbstractAnalytic::Input::Type RMT::Input::type(int index) const
{
   EDEBUG_FUNC(this,index);

   switch (index)
   {
   case InputData: return Type::DataIn;
   case LogFile: return Type::FileOut;
   case ThresholdStart: return Type::Double;
   case ThresholdStep: return Type::Double;
   case ThresholdStop: return Type::Double;
   case SplineInterpolation: return Type::Boolean;
   case MinSplinePace: return Type::Integer;
   case MaxSplinePace: return Type::Integer;
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
   EDEBUG_FUNC(this,index,role);

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
   case SplineInterpolation:
      switch (role)
      {
      case Role::CommandLineName: return QString("spline");
      case Role::Title: return tr("Use Spline Interpolation:");
      case Role::WhatsThis: return tr("Whether to perform spline interpolation on each set of eigenvalues.");
      case Role::Default: return true;
      default: return QVariant();
      }
   case MinSplinePace:
      switch (role)
      {
      case Role::CommandLineName: return QString("minpace");
      case Role::Title: return tr("Minimum Spline Pace:");
      case Role::WhatsThis: return tr("The minimum pace of the spline interpolation.");
      case Role::Default: return 10;
      case Role::Minimum: return 1;
      case Role::Maximum: return std::numeric_limits<int>::max();
      default: return QVariant();
      }
   case MaxSplinePace:
      switch (role)
      {
      case Role::CommandLineName: return QString("maxpace");
      case Role::Title: return tr("Maximum Spline Pace:");
      case Role::WhatsThis: return tr("The maximum pace of the spline interpolation.");
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
   EDEBUG_FUNC(this,index,value);

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
   case SplineInterpolation:
      _base->_splineInterpolation = value.toBool();
      break;
   case MinSplinePace:
      _base->_minSplinePace = value.toInt();
      break;
   case MaxSplinePace:
      _base->_maxSplinePace = value.toInt();
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
   EDEBUG_FUNC(this,index,file);

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
   EDEBUG_FUNC(this,index,data);

   if ( index == InputData )
   {
      _base->_input = data->cast<CorrelationMatrix>();
   }
}
