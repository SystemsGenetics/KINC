#ifndef ANALYTICFACTORY_H
#define ANALYTICFACTORY_H
#include <ace/core/AceCore.h>



class AnalyticFactory : public EAbstractAnalyticFactory
{
public:
   enum Types
   {
      ImportExpressionMatrixType = 0
      ,SpearmanType
      ,PearsonType
      ,RMTType
      ,Total
   };
   virtual quint16 getCount() override final { return Total; }
   virtual QString getName(quint16 type) override final;
   virtual QString getCommandName(quint16 type) override final;
   std::unique_ptr<EAbstractAnalytic> make(quint16 type) override final;
};



#endif
