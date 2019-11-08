#include "correlationmatrix_model.h"
#include "correlationmatrix_pair.h"



using namespace std;






/*!
 * Construct a table model for a correlation matrix.
 *
 * @param matrix
 */
CorrelationMatrix::Model::Model(CorrelationMatrix* matrix):
   _matrix(matrix)
{
   EDEBUG_FUNC(this,matrix);

   setParent(matrix);
}






/*!
 * Return a header name for the table model using a given index.
 *
 * @param section
 * @param orientation
 * @param role
 */
QVariant CorrelationMatrix::Model::headerData(int section, Qt::Orientation orientation, int role) const
{
   EDEBUG_FUNC(this,section,orientation,role);

   // orientation is not used
   Q_UNUSED(orientation);

   // if role is not display return nothing
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // get gene names
   EMetaArray geneNames {_matrix->geneNames()};

   // make sure section is within limits of gene name array
   if ( section >= 0 && section < geneNames.size() )
   {
     // return gene name
     return geneNames.at(section).toString();
   }

   // no gene found return nothing
   return QVariant();
}






/*!
 * Return the number of rows in the table model.
 *
 * @param index
 */
int CorrelationMatrix::Model::rowCount(const QModelIndex&) const
{
   EDEBUG_FUNC(this);

   return _matrix->geneSize();
}






/*!
 * Return the number of columns in the table model.
 *
 * @param index
 */
int CorrelationMatrix::Model::columnCount(const QModelIndex&) const
{
   EDEBUG_FUNC(this);

   return _matrix->geneSize();
}






/*!
 * Return a data element in the table model using the given index.
 *
 * @param index
 * @param role
 */
QVariant CorrelationMatrix::Model::data(const QModelIndex& index, int role) const
{
   EDEBUG_FUNC(this,&index,role);

   // if role is not display return nothing
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // if row and column are equal return empty string
   if ( index.row() == index.column() )
   {
      return "";
   }

   // if pairwise matrix is empty return empty string
   if ( _matrix->size() == 0 )
   {
       return "";
   }

   // get constant pair and read in values
   const Pair pair(_matrix);
   int x {index.row()};
   int y {index.column()};
   if ( y > x )
   {
      swap(x,y);
   }
   pair.read({x,y});

   // Return value of pair as a string
   return pair.toString();
}
