#ifndef ANALYTICFACTORY_H
#define ANALYTICFACTORY_H
#include <ace/core/core.h>



/*!
 * This class implements the ACE analytic factory for producing new analytic
 * objects and giving basic information about all available analytic types.
 */
class AnalyticFactory : public EAbstractAnalyticFactory
{
public:
   /*!
    * Defines all available analytic types this program implements along with the
    * total size.
    */
   enum Type
   {
      ImportExpressionMatrixType = 0
      ,ExportExpressionMatrixType
      ,ImportCorrelationMatrixType
      ,ExportCorrelationMatrixType
      ,SimilarityType
      ,CorrelationPowerFilterType
      ,ConditionalTestType
      ,PowerLawType
      ,RMTType
      ,ExtractType
      ,Total
   };
   virtual quint16 size() const override final;
   virtual QString name(quint16 type) const override final;
   virtual QString commandName(quint16 type) const override final;
   virtual std::unique_ptr<EAbstractAnalytic> make(quint16 type) const override final;
};



#endif
