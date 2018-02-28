#include "analyticfactory.h"
#include "importexpressionmatrix.h"
#include "exportexpressionmatrix.h"
#include "importcorrelationmatrix.h"
#include "exportcorrelationmatrix.h"
#include "kmeans.h"
#include "gmm.h"
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
   case ExportExpressionMatrixType: return QObject::tr("Export Expression Matrix");
   case ImportCorrelationMatrixType: return QObject::tr("Import Correlation Matrix");
   case ExportCorrelationMatrixType: return QObject::tr("Export Correlation Matrix");
   case KMeansType: return QObject::tr("K-means Clustering");
   case GMMType: return QObject::tr("GMM Clustering");
   case SpearmanType: return QObject::tr("Spearman Correlation");
   case PearsonType: return QObject::tr("Pearson Correlation");
   case RMTType: return QObject::tr("RMT Thresholding");
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
   case KMeansType: return QString("kmeans");
   case GMMType: return QString("gmm");
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
   case ExportExpressionMatrixType: return unique_ptr<EAbstractAnalytic>(new ExportExpressionMatrix);
   case ImportCorrelationMatrixType: return unique_ptr<EAbstractAnalytic>(new ImportCorrelationMatrix);
   case ExportCorrelationMatrixType: return unique_ptr<EAbstractAnalytic>(new ExportCorrelationMatrix);
   case KMeansType: return unique_ptr<EAbstractAnalytic>(new KMeans);
   case GMMType: return unique_ptr<EAbstractAnalytic>(new GMM);
   case SpearmanType: return unique_ptr<EAbstractAnalytic>(new Spearman);
   case PearsonType: return unique_ptr<EAbstractAnalytic>(new Pearson);
   case RMTType: return unique_ptr<EAbstractAnalytic>(new RMT);
   default: return nullptr;
   }
}
