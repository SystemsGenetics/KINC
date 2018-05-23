#include "datafactory.h"
#include "expressionmatrix.h"
#include "ccmatrix.h"
#include "correlationmatrix.h"



using namespace std;






quint16 DataFactory::size() const
{
   return Total;
}






QString DataFactory::name(quint16 type) const
{
   switch (type)
   {
   case ExpressionMatrixType: return "Expression Matrix";
   case CCMatrixType: return "Cluster Matrix";
   case CorrelationMatrixType: return "Correlation Matrix";
   default: return QString();
   }
}






QString DataFactory::fileExtension(quint16 type) const
{
   switch (type)
   {
   case ExpressionMatrixType: return "emx";
   case CCMatrixType: return "ccm";
   case CorrelationMatrixType: return "cmx";
   default: return QString();
   }
}






unique_ptr<EAbstractData> DataFactory::make(quint16 type) const
{
   switch (type)
   {
   case ExpressionMatrixType: return unique_ptr<EAbstractData>(new ExpressionMatrix);
   case CCMatrixType: return unique_ptr<EAbstractData>(new CCMatrix);
   case CorrelationMatrixType: return unique_ptr<EAbstractData>(new CorrelationMatrix);
   default: return nullptr;
   }
}
