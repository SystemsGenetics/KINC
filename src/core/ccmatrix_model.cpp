#include "ccmatrix_model.h"
#include "ccmatrix_pair.h"



using namespace std;






CCMatrix::Model::Model(CCMatrix* matrix):
   _matrix(matrix)
{
   setParent(matrix);
}






QVariant CCMatrix::Model::headerData(int section, Qt::Orientation orientation, int role) const
{
   // orientation is not used
   Q_UNUSED(orientation);

   // if role is not display return nothing
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // get genes metadata and make sure it is an array
   const EMetadata& genes {_matrix->geneNames()};
   if ( genes.isArray() )
   {
      // make sure section is within limits of gene name array
      if ( section >= 0 && section < genes.toArray().size() )
      {
         // return gene name
         return genes.toArray().at(section).toString();
      }
   }

   // no gene found return nothing
   return QVariant();
}






int CCMatrix::Model::rowCount(const QModelIndex&) const
{
   return _matrix->geneSize();
}






int CCMatrix::Model::columnCount(const QModelIndex&) const
{
   return _matrix->geneSize();
}






QVariant CCMatrix::Model::data(const QModelIndex &index, int role) const
{
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
