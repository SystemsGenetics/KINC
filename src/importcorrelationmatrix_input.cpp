#include "importcorrelationmatrix_input.h"
#include "datafactory.h"
#include "pairwise_index.h"



ImportCorrelationMatrix::Input::Input(ImportCorrelationMatrix* parent):
   EAbstractAnalytic::Input(parent),
   _base(parent)
{}






int ImportCorrelationMatrix::Input::size() const
{
   return Total;
}






EAbstractAnalytic::Input::Type ImportCorrelationMatrix::Input::type(int index) const
{
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






QVariant ImportCorrelationMatrix::Input::data(int index, Role role) const
{
   switch (index)
   {
   case InputFile:
      switch (role)
      {
      case Role::CommandLineName: return QString("input");
      case Role::Title: return tr("Input File:");
      case Role::WhatsThis: tr("Input text file containing pairwise correlation data.");
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
      case Role::CommandLineName: QString("genes");
      case Role::Title: return tr("Gene Size:");
      case Role::WhatsThis: return tr("Number of genes.");
      case Role::Minimum: return 1;
      case Role::Maximum: return std::numeric_limits<int>::max();
      default: return QVariant();
      }
   case MaxClusterSize:
      switch (role)
      {
      case Role::CommandLineName: QString("maxclusters");
      case Role::Title: return tr("Maximum Cluster Size:");
      case Role::WhatsThis: return tr("Maximum number of clusters per pair.");
      case Role::Minimum: return 1;
      case Role::Maximum: return Pairwise::Index::MAX_CLUSTER_SIZE;
      default: return QVariant();
      }
   case SampleSize:
      switch (role)
      {
      case Role::CommandLineName: QString("samples");
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






void ImportCorrelationMatrix::Input::set(int index, const QVariant& value)
{
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






void ImportCorrelationMatrix::Input::set(int index, QFile* file)
{
   if ( index == InputFile )
   {
      _base->_input = file;
   }
}






void ImportCorrelationMatrix::Input::set(int index, EAbstractData* data)
{
   if ( index == ClusterData )
   {
      _base->_ccm = qobject_cast<CCMatrix*>(data);
   }
   else if ( index == CorrelationData )
   {
      _base->_cmx = qobject_cast<CorrelationMatrix*>(data);
   }
}
