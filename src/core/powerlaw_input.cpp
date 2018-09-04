#include "powerlaw_input.h"
#include "correlationmatrix.h"
#include "datafactory.h"






PowerLaw::Input::Input(PowerLaw* parent):
   EAbstractAnalytic::Input(parent),
   _base(parent)
{}






int PowerLaw::Input::size() const
{
   return Total;
}






EAbstractAnalytic::Input::Type PowerLaw::Input::type(int index) const
{
   switch (index)
   {
   case InputData: return Type::DataIn;
   case LogFile: return Type::FileOut;
   case ThresholdStart: return Type::Double;
   case ThresholdStep: return Type::Double;
   case ThresholdStop: return Type::Double;
   default: return Type::Boolean;
   }
}






QVariant PowerLaw::Input::data(int index, Role role) const
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
      case Role::Default: return 0.01;
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
   default: return QVariant();
   }
}






void PowerLaw::Input::set(int index, const QVariant& value)
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
   }
}






void PowerLaw::Input::set(int index, QFile* file)
{
   if ( index == LogFile )
   {
      _base->_logfile = file;
   }
}






void PowerLaw::Input::set(int index, EAbstractData* data)
{
   if ( index == InputData )
   {
      _base->_input = data->cast<CorrelationMatrix>();
   }
}
