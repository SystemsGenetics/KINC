#include "analyticfactory.h"
#include "importexpressionmatrix.h"



using namespace std;






QString AnalyticFactory::getName(quint16 type)
{
   // figure out what data type is being queried and return name
   switch (type)
   {
   case ImportExpressionMatrixType: return QObject::tr("Import Expression Matrix");
   default: return QString();
   }
}






QString AnalyticFactory::getCommandName(quint16 type)
{
   // figure out what data type is being queried and return command name
   switch (type)
   {
   case ImportExpressionMatrixType: return QString("import_emx");
   default: return QString();
   }
}






std::unique_ptr<EAbstractAnalytic> AnalyticFactory::make(quint16 type)
{
   // figure out which data type is being requested and return new object
   switch (type)
   {
   case ImportExpressionMatrixType: return unique_ptr<EAbstractAnalytic>(new ImportExpressionMatrix);
   default: return nullptr;
   }
}
