#include "exportexpressionmatrix_input.h"
#include "datafactory.h"






ExportExpressionMatrix::Input::Input(ExportExpressionMatrix* parent):
   EAbstractAnalytic::Input(parent),
   _base(parent)
{}






int ExportExpressionMatrix::Input::size() const
{
   return Total;
}






EAbstractAnalytic::Input::Type ExportExpressionMatrix::Input::type(int index) const
{
   switch (index)
   {
   case InputData: return Type::DataIn;
   case OutputFile: return Type::FileOut;
   case NoSampleToken: return Type::String;
   default: return Type::Boolean;
   }
}






QVariant ExportExpressionMatrix::Input::data(int index, Role role) const
{
   switch (index)
   {
   case InputData:
      switch (role)
      {
      case Role::CommandLineName: return QString("input");
      case Role::Title: return tr("Input:");
      case Role::WhatsThis: return tr("Input expression matrix containing expression data.");
      case Role::DataType: return DataFactory::ExpressionMatrixType;
      default: return QVariant();
      }
   case OutputFile:
      switch (role)
      {
      case Role::CommandLineName: return QString("output");
      case Role::Title: return tr("Output:");
      case Role::WhatsThis: return tr("Output text file that will contain space/tab divided gene expression data.");
      case Role::FileFilters: return tr("Text file %1").arg("(*.txt)");
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
   default: return QVariant();
   }
}






void ExportExpressionMatrix::Input::set(int index, const QVariant& value)
{
   switch (index)
   {
   case NoSampleToken:
      _base->_noSampleToken = value.toString();
      break;
   }
}






void ExportExpressionMatrix::Input::set(int index, EAbstractData* data)
{
   if ( index == InputData )
   {
      _base->_input = data->cast<ExpressionMatrix>();
   }
}






void ExportExpressionMatrix::Input::set(int index, QFile* file)
{
   if ( index == OutputFile )
   {
      _base->_output = file;
   }
}
