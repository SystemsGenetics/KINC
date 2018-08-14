#include "extract_input.h"
#include "datafactory.h"



Extract::Input::Input(Extract* parent):
   EAbstractAnalytic::Input(parent),
   _base(parent)
{}






int Extract::Input::size() const
{
   return Total;
}






EAbstractAnalytic::Input::Type Extract::Input::type(int index) const
{
   switch (index)
   {
   case ExpressionData: return Type::DataIn;
   case ClusterData: return Type::DataIn;
   case CorrelationData: return Type::DataIn;
   case OutputFile: return Type::FileOut;
   case GraphMLFile: return Type::FileOut;
   case MinCorrelation: return Type::Double;
   case MaxCorrelation: return Type::Double;
   default: return Type::Boolean;
   }
}






QVariant Extract::Input::data(int index, Role role) const
{
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
   case CorrelationData:
      switch (role)
      {
      case Role::CommandLineName: return QString("cmx");
      case Role::Title: return tr("Correlation Matrix:");
      case Role::WhatsThis: return tr("Input correlation matrix containing correlation data.");
      case Role::DataType: return DataFactory::CorrelationMatrixType;
      default: return QVariant();
      }
   case OutputFile:
      switch (role)
      {
      case Role::CommandLineName: return QString("output");
      case Role::Title: return tr("Output File:");
      case Role::WhatsThis: return tr("Output text file that will contain network edges.");
      case Role::FileFilters: return tr("Text file %1").arg("(*.txt)");
      default: return QVariant();
      }
   case GraphMLFile:
      switch (role)
      {
      case Role::CommandLineName: return QString("graphml");
      case Role::Title: return tr("GraphML File:");
      case Role::WhatsThis: return tr("Output text file that will contain network in GraphML format.");
      case Role::FileFilters: return tr("GraphML file %1").arg("(*.graphml)");
      default: return QVariant();
      }
   case MinCorrelation:
      switch (role)
      {
      case Role::CommandLineName: return QString("mincorr");
      case Role::Title: return tr("Minimum Correlation:");
      case Role::WhatsThis: return tr("Minimum (absolute) correlation threshold for gene pairs.");
      case Role::Default: return 0.85;
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
   default: return QVariant();
   }
}






void Extract::Input::set(int index, const QVariant& value)
{
   switch (index)
   {
   case MinCorrelation:
      _base->_minCorrelation = value.toFloat();
      break;
   case MaxCorrelation:
      _base->_maxCorrelation = value.toFloat();
      break;
   }
}






void Extract::Input::set(int index, EAbstractData* data)
{
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






void Extract::Input::set(int index, QFile* file)
{
   if ( index == OutputFile )
   {
      _base->_output = file;
   }
   else if ( index == GraphMLFile )
   {
      _base->_graphml = file;
   }
}
