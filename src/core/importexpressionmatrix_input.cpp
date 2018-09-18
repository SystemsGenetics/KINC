#include "importexpressionmatrix_input.h"
#include "datafactory.h"






ImportExpressionMatrix::Input::Input(ImportExpressionMatrix* parent):
   EAbstractAnalytic::Input(parent),
   _base(parent)
{}






int ImportExpressionMatrix::Input::size() const
{
   return Total;
}






EAbstractAnalytic::Input::Type ImportExpressionMatrix::Input::type(int index) const
{
   switch (index)
   {
   case InputFile: return Type::FileIn;
   case OutputData: return Type::DataOut;
   case NoSampleToken: return Type::String;
   case SampleSize: return Type::Integer;
   default: return Type::Boolean;
   }
}






QVariant ImportExpressionMatrix::Input::data(int index, Role role) const
{
   switch (index)
   {
   case InputFile:
      switch (role)
      {
      case Role::CommandLineName: return QString("input");
      case Role::Title: return tr("Input:");
      case Role::WhatsThis: return tr("Input text file containing space/tab delimited gene expression data.");
      case Role::FileFilters: return tr("Text file %1").arg("(*.txt)");
      default: return QVariant();
      }
   case OutputData:
      switch (role)
      {
      case Role::CommandLineName: return QString("output");
      case Role::Title: return tr("Output:");
      case Role::WhatsThis: return tr("Output expression matrix that will contain expression data.");
      case Role::DataType: return DataFactory::ExpressionMatrixType;
      default: return QVariant();
      }
   case NoSampleToken:
      switch (role)
      {
      case Role::CommandLineName: return QString("nan");
      case Role::Title: return tr("No Sample Token:");
      case Role::WhatsThis: return tr("Expected token for expressions that have no value.");
      default: return QVariant();
      }
   case SampleSize:
      switch (role)
      {
      case Role::CommandLineName: return QString("samples");
      case Role::Title: return tr("Sample Size:");
      case Role::WhatsThis: return tr("Number of samples. 0 indicates the text file contains a header of sample names to be read to determine size.");
      case Role::Default: return 0;
      case Role::Minimum: return 0;
      case Role::Maximum: return std::numeric_limits<int>::max();
      default: return QVariant();
      }
   default: return QVariant();
   }
}






void ImportExpressionMatrix::Input::set(int index, const QVariant& value)
{
   switch (index)
   {
   case SampleSize:
      _base->_sampleSize = value.toInt();
      break;
   case NoSampleToken:
      _base->_noSampleToken = value.toString();
      break;
   }
}






void ImportExpressionMatrix::Input::set(int index, QFile* file)
{
   if ( index == InputFile )
   {
      _base->_input = file;
   }
}






void ImportExpressionMatrix::Input::set(int index, EAbstractData* data)
{
   if ( index == OutputData )
   {
      _base->_output = data->cast<ExpressionMatrix>();
   }
}
