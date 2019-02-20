#include "expressionmatrix_model.h"
//






/*!
 * Construct a table model for an expression matrix.
 *
 * @param matrix
 */
ExpressionMatrix::Model::Model(ExpressionMatrix* matrix):
   _matrix(matrix)
{
   EDEBUG_FUNC(this,matrix);

   setParent(matrix);
}






/*!
 * Return a header name for the table model using a given index and
 * orientation (row / column).
 *
 * @param section
 * @param orientation
 * @param role
 */
QVariant ExpressionMatrix::Model::headerData(int section, Qt::Orientation orientation, int role) const
{
   EDEBUG_FUNC(this,section,orientation,role);

   // make sure the role is valid
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // determine whether to return a row name or column name
   switch (orientation)
   {
   case Qt::Vertical:
   {
      // get gene names
      EMetaArray geneNames {_matrix->geneNames()};

      // make sure the index is valid
      if ( section >= 0 && section < geneNames.size() )
      {
         // return the specified row name
         return geneNames.at(section).toString();
      }

      // otherwise return empty string
      return QVariant();
   }
   case Qt::Horizontal:
   {
      // get sample names
      EMetaArray samples {_matrix->sampleNames()};

      // make sure the index is valid
      if ( section >= 0 && section < samples.size() )
      {
         // return the specified column name
         return samples.at(section).toString();
      }

      // otherwise return empty string
      return QVariant();
   }
   }

   return QVariant();
}






/*!
 * Return the number of rows in the table model.
 *
 * @param index
 */
int ExpressionMatrix::Model::rowCount(const QModelIndex&) const
{
   EDEBUG_FUNC(this);

   return _matrix->_geneSize;
}






/*!
 * Return the number of columns in the table model.
 *
 * @param index
 */
int ExpressionMatrix::Model::columnCount(const QModelIndex&) const
{
   EDEBUG_FUNC(this);

   return _matrix->_sampleSize;
}






/*!
 * Return a data element in the table model using the given index.
 *
 * @param index
 * @param role
 */
QVariant ExpressionMatrix::Model::data(const QModelIndex& index, int role) const
{
   EDEBUG_FUNC(this,&index,role);

   // make sure the index and role are valid
   if ( !index.isValid() || role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // make sure the index is within the bounds of the expression matrix
   if ( index.row() >= _matrix->_geneSize || index.column() >= _matrix->_sampleSize )
   {
      return QVariant();
   }

   // get the specified value from the expression matrix
   float value;
   _matrix->seekExpression(index.row(),index.column());
   _matrix->stream() >> value;

   // return the specified value
   return value;
}
