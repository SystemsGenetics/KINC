#include "datafactory.h"
#include "expressionmatrix.h"



using namespace std;






QString DataFactory::getName(quint16 type) noexcept
{
   // figure out what data type is being queried and return name
   switch (type)
   {
   case ExpressionMatrixType: return QObject::tr("Expression Matrix");
   default: return QString();
   }
}






QString DataFactory::getFileExtension(quint16 type) noexcept
{
   // figure out what data type is being queried and return extension
   switch (type)
   {
   case ExpressionMatrixType: return QString("emx");
   default: return QString();
   }
}






unique_ptr<EAbstractData> DataFactory::make(quint16 type) noexcept
{
   // figure out which data type is being requested and return new object
   switch (type)
   {
   case ExpressionMatrixType: return unique_ptr<EAbstractData>(new ExpressionMatrix);
   default: return nullptr;
   }
}
