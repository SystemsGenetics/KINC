#ifndef DATAFACTORY_H
#define DATAFACTORY_H
#include <ace/core/core.h>



class DataFactory : public EAbstractDataFactory
{
public:
   enum Type
   {
      ExpressionMatrixType = 0
      ,CCMatrixType
      ,CorrelationMatrixType
      ,Total
   };
   virtual quint16 size() const override final;
   virtual QString name(quint16 type) const override final;
   virtual QString fileExtension(quint16 type) const override final;
   virtual std::unique_ptr<EAbstractData> make(quint16 type) const override final;
};



#endif
