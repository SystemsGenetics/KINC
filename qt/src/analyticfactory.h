#ifndef ANALYTICFACTORY_H
#define ANALYTICFACTORY_H
#include <ace/core/AceCore.h>



class AnalyticFactory : public EAbstractAnalyticFactory
{
public:
   enum Types
   {
      Total = 0
   };
   virtual quint16 getCount() override final { return Total; }
   virtual QString getName(quint16 type) override final { return QString(); }
   virtual QString getCommandName(quint16 type) override final { return QString(); }
   std::unique_ptr<EAbstractAnalytic> make(quint16 type) { return nullptr; }
};



#endif
