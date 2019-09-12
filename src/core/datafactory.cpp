#include "datafactory.h"
#include "expressionmatrix.h"
#include "ccmatrix.h"
#include "correlationmatrix.h"
#include "conditionspecificclustersmatrix.h"



using namespace std;






/*!
 * Return the total number of data types this program implements.
 */
quint16 DataFactory::size() const
{
   EDEBUG_FUNC(this);

   return Total;
}






/*!
 * Return the display name for the given data type.
 *
 * @param type
 */
QString DataFactory::name(quint16 type) const
{
   EDEBUG_FUNC(this,type);

   switch (type)
   {
   case ExpressionMatrixType: return "Expression Matrix";
   case CCMatrixType: return "Cluster Matrix";
   case CorrelationMatrixType: return "Correlation Matrix";
   case CSMType: return "Condition-Specific Clusters Matrix";
   default: return QString();
   }
}






/*!
 * Return the file extension for the given data type as a string.
 *
 * @param type
 */
QString DataFactory::fileExtension(quint16 type) const
{
   EDEBUG_FUNC(this,type);

   switch (type)
   {
   case ExpressionMatrixType: return "emx";
   case CCMatrixType: return "ccm";
   case CorrelationMatrixType: return "cmx";
   case CSMType: return "csm";
   default: return QString();
   }
}






/*!
 * Make and return a new abstract data object of the given type.
 *
 * @param type
 */
unique_ptr<EAbstractData> DataFactory::make(quint16 type) const
{
   EDEBUG_FUNC(this,type);
   switch (type)
   {
   case ExpressionMatrixType: return unique_ptr<EAbstractData>(new ExpressionMatrix);
   case CCMatrixType: return unique_ptr<EAbstractData>(new CCMatrix);
   case CorrelationMatrixType: return unique_ptr<EAbstractData>(new CorrelationMatrix);
   case CSMType: return unique_ptr<CSM>(new CSM);
   default: return nullptr;
   }
}
