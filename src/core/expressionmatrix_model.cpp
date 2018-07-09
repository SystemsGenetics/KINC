#include "expressionmatrix_model.h"
//






/*!
 *
 * @param matrix  
 */
ExpressionMatrix::Model::Model(ExpressionMatrix* matrix):
   _matrix(matrix)
{
   setParent(matrix);
}






/*!
 *
 * @param section  
 *
 * @param orientation  
 *
 * @param role  
 */
QVariant ExpressionMatrix::Model::headerData(int section, Qt::Orientation orientation, int role) const
{
   // if this is not display role return nothing
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // get metadata root and figure out orientation
   switch (orientation)
   {
   case Qt::Vertical:
   {
      // get gene names and make sure it is array
      EMetadata genes {_matrix->meta().toObject().at("genes")};
      if ( genes.isArray() )
      {
         // make sure section is within limits of array
         if ( section >= 0 && section < genes.toArray().size() )
         {
            // return gene name
            return genes.toArray().at(section).toString();
         }
      }

      // if no gene name found return nothing
      return QVariant();
   }
   case Qt::Horizontal:
   {
      // get sample names and make sure it is array
      EMetadata samples {_matrix->meta().toObject().at("samples")};
      if ( samples.isArray() )
      {
         // make sure section is within limits of array
         if ( section >= 0 && section < samples.toArray().size() )
         {
            // return sample name
            return samples.toArray().at(section).toString();
         }
      }

      // if no sample name found return nothing
      return QVariant();
   }
   default:
      // unknown orientation so return nothing
      return QVariant();
   }
}






/*!
 *
 * @param parent  
 */
int ExpressionMatrix::Model::rowCount(const QModelIndex& parent) const
{
   // return gene size for row count
   Q_UNUSED(parent);
   return _matrix->_geneSize;
}






/*!
 *
 * @param parent  
 */
int ExpressionMatrix::Model::columnCount(const QModelIndex& parent) const
{
   // return sample size for column count
   Q_UNUSED(parent);
   return _matrix->_sampleSize;
}






/*!
 *
 * @param index  
 *
 * @param role  
 */
QVariant ExpressionMatrix::Model::data(const QModelIndex& index, int role) const
{
   // if role is not display return nothing
   if ( !index.isValid() || role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // if index is out of range return nothing
   if ( index.row() >= _matrix->_geneSize || index.column() >= _matrix->_sampleSize )
   {
      return QVariant();
   }

   // make input variable and seek to position of queried expression
   float value;
   _matrix->seekExpression(index.row(),index.column());
   _matrix->stream() >> value;

   // return expression
   return value;
}
