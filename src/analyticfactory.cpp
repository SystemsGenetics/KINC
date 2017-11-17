#include "analyticfactory.h"
#include "importexpressionmatrix.h"
#include "spearman.h"
#include "pearson.h"
#include "rmt.h"



using namespace std;






QString AnalyticFactory::getName(quint16 type)
{
   // figure out what data type is being queried and return name
   switch (type)
   {
   case ImportExpressionMatrixType: return QObject::tr("Import Expression Matrix");
   case SpearmanType: return QObject::tr("Spearman");
   case PearsonType: return QObject::tr("Pearson");
   case RMTType: return QObject::tr("RMT");
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
   case PearsonType: return QString("pearson");
   case RMTType: return QString("rmt");
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
   case PearsonType: return unique_ptr<EAbstractAnalytic>(new Pearson);
   case RMTType: return unique_ptr<EAbstractAnalytic>(new RMT);
   default: return nullptr;
   }
}
