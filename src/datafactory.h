#ifndef DATAFACTORY_H
#define DATAFACTORY_H
#include <ace/core/AceCore.h>



class DataFactory : public EAbstractDataFactory
{
public:
   enum Types
   {
      ExpressionMatrixType = 0
      ,CorrelationMatrixType
      ,Total
   };
   virtual quint16 getCount() override final { return Total; }
   virtual QString getName(quint16 type) override final;
   virtual QString getFileExtension(quint16 type) override final;
   virtual std::unique_ptr<EAbstractData> make(quint16 type) override final;
};



#endif
