#include "analyticfactory.h"
#include "importexpressionmatrix.h"
#include "exportexpressionmatrix.h"
#include "importcorrelationmatrix.h"
#include "exportcorrelationmatrix.h"
#include "similarity.h"
#include "rmt.h"
#include "extract.h"



using namespace std;






QString AnalyticFactory::getName(quint16 type)
{
   // figure out what data type is being queried and return name
   switch (type)
   {
   case ImportExpressionMatrixType: return QObject::tr("Import Expression Matrix");
   case ExportExpressionMatrixType: return QObject::tr("Export Expression Matrix");
   case ImportCorrelationMatrixType: return QObject::tr("Import Correlation Matrix");
   case ExportCorrelationMatrixType: return QObject::tr("Export Correlation Matrix");
   case SimilarityType: return QObject::tr("Similarity");
   case RMTType: return QObject::tr("RMT Thresholding");
   case ExtractType: return QObject::tr("Extract Network");
   default: return QString();
   }
}






QString AnalyticFactory::getCommandName(quint16 type)
{
   // figure out what data type is being queried and return command name
   switch (type)
   {
   case ImportExpressionMatrixType: return QString("import_emx");
   case ExportExpressionMatrixType: return QString("export_emx");
   case ImportCorrelationMatrixType: return QString("import_cmx");
   case ExportCorrelationMatrixType: return QString("export_cmx");
   case SimilarityType: return QString("similarity");
   case RMTType: return QString("rmt");
   case ExtractType: return QString("extract");
   default: return QString();
   }
}






std::unique_ptr<EAbstractAnalytic> AnalyticFactory::make(quint16 type)
{
   // figure out which data type is being requested and return new object
   switch (type)
   {
   case ImportExpressionMatrixType: return unique_ptr<EAbstractAnalytic>(new ImportExpressionMatrix);
   case ExportExpressionMatrixType: return unique_ptr<EAbstractAnalytic>(new ExportExpressionMatrix);
   case ImportCorrelationMatrixType: return unique_ptr<EAbstractAnalytic>(new ImportCorrelationMatrix);
   case ExportCorrelationMatrixType: return unique_ptr<EAbstractAnalytic>(new ExportCorrelationMatrix);
   case SimilarityType: return unique_ptr<EAbstractAnalytic>(new Similarity);
   case RMTType: return unique_ptr<EAbstractAnalytic>(new RMT);
   case ExtractType: return unique_ptr<EAbstractAnalytic>(new Extract);
   default: return nullptr;
   }
}
