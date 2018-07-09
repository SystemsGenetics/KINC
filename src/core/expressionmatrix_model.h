#ifndef EXPRESSIONMATRIX_MODEL_H
#define EXPRESSIONMATRIX_MODEL_H
#include "expressionmatrix.h"
//



/*!
 */
class ExpressionMatrix::Model : public QAbstractTableModel
{
public:
   Model(ExpressionMatrix* matrix);
   virtual QVariant headerData(int section, Qt::Orientation orientation, int role) const override final;
   virtual int rowCount(const QModelIndex& parent) const override final;
   virtual int columnCount(const QModelIndex& parent) const override final;
   virtual QVariant data(const QModelIndex& index, int role) const override final;
private:
   /*!
    */
   ExpressionMatrix* _matrix;
};



#endif
