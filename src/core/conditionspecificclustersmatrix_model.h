#ifndef CSM_MODEL_H
#define CSM_MODEL_H
#include "conditionspecificclustersmatrix.h"
#include <ace/core/core.h>



class CSMatrix::Model : public QAbstractTableModel
{
public:
    Model(CSMatrix* matrix);
public:
    virtual QVariant headerData(int section, Qt::Orientation orientation, int role) const override final;
    virtual int columnCount(const QModelIndex&) const override final;
    virtual int rowCount(const QModelIndex&) const override final;
    virtual QVariant data(const QModelIndex& index, int role) const override final;
private:
    /*!
     * Pointer to the data object for this table model.
     */
    CSMatrix* _matrix;
};



#endif
