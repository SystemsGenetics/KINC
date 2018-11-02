#include "pairwise_index.h"




using namespace Pairwise;






/*!
 * Construct a pairwise index from a row index and a column index. The row
 * index must be greater than the column index.
 *
 * @param x
 * @param y
 */
Index::Index(qint32 x, qint32 y):
   _x(x),
   _y(y)
{
   EDEBUG_FUNC(this,x,y);

   // make sure pairwise index is valid
   if ( x < 1 || y < 0 || x <= y )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Pairwise Index Error"));
      e.setDetails(QObject::tr("Cannot initialize pairwise index (%1, %2).").arg(x).arg(y));
      throw e;
   }
}






/*!
 * Construct a pairwise index from a one-dimensional index, which corresponds
 * to the i-th element in the lower triangle of a matrix using row-major order.
 *
 * @param index
 */
Index::Index(qint64 index)
{
   EDEBUG_FUNC(this,index);

   // make sure index is valid
   if ( index < 0 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Pairwise Index Error"));
      e.setDetails(QObject::tr("Cannot initialize pairwise index from %1.").arg(index));
      throw e;
   }

   // compute pairwise index from scalar index
   qint64 pos {0};
   qint64 x {0};

   while ( pos + x <= index )
   {
      pos += x;
      ++x;
   }

   _x = x;
   _y = index - pos;
}






/*!
 * Return the indent value of this pairwise index with a given cluster index.
 *
 * @param cluster
 */
qint64 Index::indent(qint8 cluster) const
{
   EDEBUG_FUNC(this,cluster);

   // make sure cluster given is valid
   if ( cluster < 0 || cluster >= MAX_CLUSTER_SIZE )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Pairwise Index Error"));
      e.setDetails(QObject::tr("Cluster %1 is outside limits.").arg(cluster));
      throw e;
   }

   // compute indent with given cluster and return it
   qint64 index {(qint64)_x * (_x - 1) / 2 + _y};
   return index * MAX_CLUSTER_SIZE + cluster;
}






/*!
 * Increment a pairwise index to the next element.
 */
void Index::operator++()
{
   EDEBUG_FUNC(this);

   // increment gene y and check if it reaches gene x
   if ( ++_y >= _x )
   {
      // reset gene y to 0 and increment gene x
      _y = 0;
      ++_x;
   }
}
