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
   // 
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // 
   switch (orientation)
   {
   case Qt::Vertical:
   {
      // 
      EMetadata genes {_matrix->meta().toObject().at("genes")};
      if ( genes.isArray() )
      {
         // 
         if ( section >= 0 && section < genes.toArray().size() )
         {
            // 
            return genes.toArray().at(section).toString();
         }
      }

      // 
      return QVariant();
   }
   case Qt::Horizontal:
   {
      // 
      EMetadata samples {_matrix->meta().toObject().at("samples")};
      if ( samples.isArray() )
      {
         // 
         if ( section >= 0 && section < samples.toArray().size() )
         {
            // 
            return samples.toArray().at(section).toString();
         }
      }

      // 
      return QVariant();
   }
   default:
      // 
      return QVariant();
   }
}






/*!
 *
 * @param parent  
 */
int ExpressionMatrix::Model::rowCount(const QModelIndex& parent) const
{
   // 
   Q_UNUSED(parent);
   return _matrix->_geneSize;
}






/*!
 *
 * @param parent  
 */
int ExpressionMatrix::Model::columnCount(const QModelIndex& parent) const
{
   // 
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
   // 
   if ( !index.isValid() || role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // 
   if ( index.row() >= _matrix->_geneSize || index.column() >= _matrix->_sampleSize )
   {
      return QVariant();
   }

   // 
   float value;
   _matrix->seekExpression(index.row(),index.column());
   _matrix->stream() >> value;

   // 
   return value;
}
