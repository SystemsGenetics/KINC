#ifndef CSCM_MODEL_H
#define CSCM_MODEL_H
#include "condition-specificclustersmatrix.h"
#include <ace/core/core.h>

class CSCM::Model : public QAbstractTableModel
{
public:
   Model(CSCM* matrix);
public:
   virtual QVariant headerData(int section, Qt::Orientation orientation, int role) const override final;
   virtual int columnCount(const QModelIndex&) const override final;
   virtual int rowCount(const QModelIndex&) const override final;
   virtual QVariant data(const QModelIndex& index, int role) const override final;
private:
   /*!
    * Pointer to the data object for this table model.
    */
   CSCM* _matrix;
};

#endif // CSCM_MODEL_H
