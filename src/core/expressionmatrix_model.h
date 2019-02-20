#ifndef EXPRESSIONMATRIX_MODEL_H
#define EXPRESSIONMATRIX_MODEL_H
#include "expressionmatrix.h"
//



/*!
 * This class implements the qt table model for the expression matrix
 * data object, which represents the expression matrix as a table.
 */
class ExpressionMatrix::Model : public QAbstractTableModel
{
public:
   Model(ExpressionMatrix* matrix);
public:
   virtual QVariant headerData(int section, Qt::Orientation orientation, int role) const override final;
   virtual int rowCount(const QModelIndex&) const override final;
   virtual int columnCount(const QModelIndex&) const override final;
   virtual QVariant data(const QModelIndex& index, int role) const override final;
private:
   /*!
    * Pointer to the data object for this table model.
    */
   ExpressionMatrix* _matrix;
};



#endif
