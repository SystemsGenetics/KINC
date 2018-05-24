#include "analyticfactory.h"
#include "importexpressionmatrix.h"
// #include "exportexpressionmatrix.h"
// #include "importcorrelationmatrix.h"
// #include "exportcorrelationmatrix.h"
// #include "similarity.h"
// #include "rmt.h"
// #include "extract.h"



using namespace std;






quint16 AnalyticFactory::size() const
{
   return Total;
}






QString AnalyticFactory::name(quint16 type) const
{
   switch (type)
   {
   case ImportExpressionMatrixType: return "Import Expression Matrix";
   // case ExportExpressionMatrixType: return "Export Expression Matrix";
   // case ImportCorrelationMatrixType: return "Import Correlation Matrix";
   // case ExportCorrelationMatrixType: return "Export Correlation Matrix";
   // case SimilarityType: return "Similarity";
   // case RMTType: return "RMT Thresholding";
   // case ExtractType: return "Extract Network";
   default: return QString();
   }
}






QString AnalyticFactory::commandName(quint16 type) const
{
   switch (type)
   {
   case ImportExpressionMatrixType: return "import-emx";
   // case ExportExpressionMatrixType: return "export_emx";
   // case ImportCorrelationMatrixType: return "import_cmx";
   // case ExportCorrelationMatrixType: return "export_cmx";
   // case SimilarityType: return "similarity";
   // case RMTType: return "rmt";
   // case ExtractType: return "extract";
   default: return QString();
   }
}






std::unique_ptr<EAbstractAnalytic> AnalyticFactory::make(quint16 type) const
{
   switch (type)
   {
   case ImportExpressionMatrixType: return unique_ptr<EAbstractAnalytic>(new ImportExpressionMatrix);
   // case ExportExpressionMatrixType: return unique_ptr<EAbstractAnalytic>(new ExportExpressionMatrix);
   // case ImportCorrelationMatrixType: return unique_ptr<EAbstractAnalytic>(new ImportCorrelationMatrix);
   // case ExportCorrelationMatrixType: return unique_ptr<EAbstractAnalytic>(new ExportCorrelationMatrix);
   // case SimilarityType: return unique_ptr<EAbstractAnalytic>(new Similarity);
   // case RMTType: return unique_ptr<EAbstractAnalytic>(new RMT);
   // case ExtractType: return unique_ptr<EAbstractAnalytic>(new Extract);
   default: return nullptr;
   }
}
