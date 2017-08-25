#ifndef DATAFACTORY_H
#define DATAFACTORY_H
#include <ace/core/AceCore.h>



class DataFactory : public EAbstractDataFactory
{
public:
   enum Types
   {
      Total = 0
   };
   virtual quint16 getCount() noexcept override final { return Total; }
   virtual QString getName(quint16 type) noexcept override final { return QString(); }
   virtual QString getFileExtension(quint16 type) noexcept override final { return QString(); }
   virtual std::unique_ptr<EAbstractData> make(quint16 type) noexcept override final
      { return nullptr; }
};



#endif
