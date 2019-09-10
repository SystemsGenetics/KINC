#ifndef CPMATRIX_MODEL_H
#define CPMATRIX_MODEL_H
#include "cpmatrix.h"



/*!
 * This class implements the qt table model for the cluster matrix
 * data object, which represents the cluster matrix as a table.
 */
class CPMatrix::Model : public QAbstractTableModel
{
public:
   Model(CPMatrix* matrix);
public:
   virtual QVariant headerData(int section, Qt::Orientation orientation, int role) const override final;
   virtual int rowCount(const QModelIndex&) const override final;
   virtual int columnCount(const QModelIndex&) const override final;
   virtual QVariant data(const QModelIndex& index, int role) const override final;
private:
   /*!
    * Pointer to the data object for this table model.
    */
   CPMatrix* _matrix;
};



#endif
