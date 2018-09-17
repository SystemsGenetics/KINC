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
      // get the specified row name from the data object's metadata
      EMetadata genes {_matrix->meta().toObject().at("genes")};
      if ( genes.isArray() )
      {
         // make sure the index is valid
         if ( section >= 0 && section < genes.toArray().size() )
         {
            // return the specified row name
            return genes.toArray().at(section).toString();
         }
      }

      // otherwise return empty string
      return QVariant();
   }
   case Qt::Horizontal:
   {
      // get the specified column name from the data object's metadata
      EMetadata samples {_matrix->meta().toObject().at("samples")};
      if ( samples.isArray() )
      {
         // make sure the index is valid
         if ( section >= 0 && section < samples.toArray().size() )
         {
            // return the specified column name
            return samples.toArray().at(section).toString();
         }
      }

      // otherwise return empty string
      return QVariant();
   }
   default:
      // return empty string if orientation is not valid
      return QVariant();
   }
}






/*!
 * Return the number of rows in the table model.
 *
 * @param parent
 */
int ExpressionMatrix::Model::rowCount(const QModelIndex& parent) const
{
   Q_UNUSED(parent);
   return _matrix->_geneSize;
}






/*!
 * Return the number of columns in the table model.
 *
 * @param parent
 */
int ExpressionMatrix::Model::columnCount(const QModelIndex& parent) const
{
   Q_UNUSED(parent);
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
