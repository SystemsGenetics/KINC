#ifndef IMPORTCORRELATIONMATRIX_H
#define IMPORTCORRELATIONMATRIX_H
#include <ace/core/core.h>

#include "ccmatrix.h"
#include "correlationmatrix.h"



class ImportCorrelationMatrix : public EAbstractAnalytic
{
   Q_OBJECT
public:
   class Input;
   virtual int size() const override final;
   virtual void process(const EAbstractAnalytic::Block* result) override final;
   virtual EAbstractAnalytic::Input* makeInput() override final;
   virtual void initialize();
private:
   QFile* _input {nullptr};
   CCMatrix* _ccm {nullptr};
   CorrelationMatrix* _cmx {nullptr};
   qint32 _geneSize {0};
   qint32 _maxClusterSize {1};
   qint32 _sampleSize {0};
   QString _correlationName;
};



#endif
