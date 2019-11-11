#include "analyticfactory.h"
#include "importexpressionmatrix.h"
#include "exportexpressionmatrix.h"
#include "importcorrelationmatrix.h"
#include "exportcorrelationmatrix.h"
#include "conditionaltest.h"
#include "similarity.h"
#include "corrpower.h"
#include "powerlaw.h"
#include "rmt.h"
#include "extract.h"



using namespace std;






/*!
 * Return the total number of analytic types that this program implements.
 */
quint16 AnalyticFactory::size() const
{
   EDEBUG_FUNC(this);

   return Total;
}






/*!
 * Return the display name for the given analytic type.
 *
 * @param type
 */
QString AnalyticFactory::name(quint16 type) const
{
   EDEBUG_FUNC(this,type);

   switch (type)
   {
   case ImportExpressionMatrixType: return "Import Expression Matrix";
   case ExportExpressionMatrixType: return "Export Expression Matrix";
   case ImportCorrelationMatrixType: return "Import Correlation Matrix";
   case ExportCorrelationMatrixType: return "Export Correlation Matrix";
   case SimilarityType: return "Similarity";
   case CorrelationPowerFilterType: return "Filter: Correlation Power";
   case ConditionalTestType: return "Threshold: Condition-Specific";
   case PowerLawType: return "Threshold: Power-law";
   case RMTType: return "Threshold: RMT";
   case ExtractType: return "Extract Network";
   default: return QString();
   }
}






/*!
 * Return the command line name for the given analytic type.
 *
 * @param type
 */
QString AnalyticFactory::commandName(quint16 type) const
{
   EDEBUG_FUNC(this,type);

   switch (type)
   {
   case ImportExpressionMatrixType: return "import-emx";
   case ExportExpressionMatrixType: return "export-emx";
   case ImportCorrelationMatrixType: return "import-cmx";
   case ExportCorrelationMatrixType: return "export-cmx";
   case SimilarityType: return "similarity";
   case CorrelationPowerFilterType: return "corrpower";
   case ConditionalTestType: return "cond-test";
   case PowerLawType: return "powerlaw";
   case RMTType: return "rmt";
   case ExtractType: return "extract";
   default: return QString();
   }
}






/*!
 * Make and return a new abstract analytic object of the given type.
 *
 * @param type
 */
std::unique_ptr<EAbstractAnalytic> AnalyticFactory::make(quint16 type) const
{
   EDEBUG_FUNC(this,type);

   switch (type)
   {
   case ImportExpressionMatrixType: return unique_ptr<EAbstractAnalytic>(new ImportExpressionMatrix);
   case ExportExpressionMatrixType: return unique_ptr<EAbstractAnalytic>(new ExportExpressionMatrix);
   case ImportCorrelationMatrixType: return unique_ptr<EAbstractAnalytic>(new ImportCorrelationMatrix);
   case ExportCorrelationMatrixType: return unique_ptr<EAbstractAnalytic>(new ExportCorrelationMatrix);
   case SimilarityType: return unique_ptr<EAbstractAnalytic>(new Similarity);
   case CorrelationPowerFilterType: return unique_ptr<EAbstractAnalytic>(new CorrPowerFilter);
   case ConditionalTestType: return unique_ptr<EAbstractAnalytic>(new ConditionalTest);
   case PowerLawType: return unique_ptr<EAbstractAnalytic>(new PowerLaw);
   case RMTType: return unique_ptr<EAbstractAnalytic>(new RMT);
   case ExtractType: return unique_ptr<EAbstractAnalytic>(new Extract);
   default: return nullptr;
   }
}
