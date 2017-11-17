#include "datafactory.h"
#include "expressionmatrix.h"
#include "correlationmatrix.h"



using namespace std;






QString DataFactory::getName(quint16 type)
{
   // figure out what data type is being queried and return name
   switch (type)
   {
   case ExpressionMatrixType: return QObject::tr("Expression Matrix");
   case CorrelationMatrixType: return QObject::tr("Correlation Matrix");
   default: return QString();
   }
}






QString DataFactory::getFileExtension(quint16 type)
{
   // figure out what data type is being queried and return extension
   switch (type)
   {
   case ExpressionMatrixType: return QString("emx");
   case CorrelationMatrixType: return QString("cmx");
   default: return QString();
   }
}






unique_ptr<EAbstractData> DataFactory::make(quint16 type)
{
   // figure out which data type is being requested and return new object
   switch (type)
   {
   case ExpressionMatrixType: return unique_ptr<EAbstractData>(new ExpressionMatrix);
   case CorrelationMatrixType: return unique_ptr<EAbstractData>(new CorrelationMatrix);
   default: return nullptr;
   }
}
