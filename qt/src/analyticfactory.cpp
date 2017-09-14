#include "analyticfactory.h"
#include "importexpressionmatrix.h"
#include "spearman.h"



using namespace std;






QString AnalyticFactory::getName(quint16 type)
{
   // figure out what data type is being queried and return name
   switch (type)
   {
   case ImportExpressionMatrixType: return QObject::tr("Import Expression Matrix");
   case SpearmanType: return QObject::tr("Spearman");
   default: return QString();
   }
}






QString AnalyticFactory::getCommandName(quint16 type)
{
   // figure out what data type is being queried and return command name
   switch (type)
   {
   case ImportExpressionMatrixType: return QString("import_emx");
   case SpearmanType: return QString("spearman");
   default: return QString();
   }
}






std::unique_ptr<EAbstractAnalytic> AnalyticFactory::make(quint16 type)
{
   // figure out which data type is being requested and return new object
   switch (type)
   {
   case ImportExpressionMatrixType: return unique_ptr<EAbstractAnalytic>(new ImportExpressionMatrix);
   case SpearmanType: return unique_ptr<EAbstractAnalytic>(new Spearman);
   default: return nullptr;
   }
}
