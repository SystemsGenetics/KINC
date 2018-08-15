#ifndef CCMATRIX_MODEL_H
#define CCMATRIX_MODEL_H
#include "ccmatrix.h"



class CCMatrix::Model : public QAbstractTableModel
{
public:
   Model(CCMatrix* matrix);
   virtual QVariant headerData(int section, Qt::Orientation orientation, int role) const override final;
   virtual int rowCount(const QModelIndex& parent) const override final;
   virtual int columnCount(const QModelIndex& parent) const override final;
   virtual QVariant data(const QModelIndex& index, int role) const override final;
private:
   CCMatrix* _matrix;
};



#endif
