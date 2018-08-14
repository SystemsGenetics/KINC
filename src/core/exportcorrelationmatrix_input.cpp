#include "exportcorrelationmatrix_input.h"
#include "datafactory.h"



ExportCorrelationMatrix::Input::Input(ExportCorrelationMatrix* parent):
   EAbstractAnalytic::Input(parent),
   _base(parent)
{}






int ExportCorrelationMatrix::Input::size() const
{
   return Total;
}






EAbstractAnalytic::Input::Type ExportCorrelationMatrix::Input::type(int index) const
{
   switch (index)
   {
   case ExpressionData: return Type::DataIn;
   case ClusterData: return Type::DataIn;
   case CorrelationData: return Type::DataIn;
   case OutputFile: return Type::FileOut;
   default: return Type::Boolean;
   }
}






QVariant ExportCorrelationMatrix::Input::data(int index, Role role) const
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
      case Role::WhatsThis: return tr("Output text file that will contain pairwise correlation data.");
      case Role::FileFilters: return tr("Text file %1").arg("(*.txt)");
      default: return QVariant();
      }
   default: return QVariant();
   }
}






void ExportCorrelationMatrix::Input::set(int index, const QVariant& value)
{
   Q_UNUSED(index);
   Q_UNUSED(value);
}






void ExportCorrelationMatrix::Input::set(int index, EAbstractData* data)
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






void ExportCorrelationMatrix::Input::set(int index, QFile* file)
{
   if ( index == OutputFile )
   {
      _base->_output = file;
   }
}
