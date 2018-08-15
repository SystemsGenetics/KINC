#ifndef CORRELATIONMATRIX_MODEL_H
#define CORRELATIONMATRIX_MODEL_H
#include "correlationmatrix.h"



class CorrelationMatrix::Model : public QAbstractTableModel
{
public:
   Model(CorrelationMatrix* matrix);
   virtual QVariant headerData(int section, Qt::Orientation orientation, int role) const override final;
   virtual int rowCount(const QModelIndex& parent) const override final;
   virtual int columnCount(const QModelIndex& parent) const override final;
   virtual QVariant data(const QModelIndex& index, int role) const override final;
private:
   CorrelationMatrix* _matrix;
};



#endif
