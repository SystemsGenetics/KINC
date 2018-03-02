#include "genepair_vector.h"




using namespace GenePair;






Vector::Vector(qint32 geneX, qint32 geneY):
   _geneX(geneX),
   _geneY(geneY)
{
   // make sure vector given is valid
   if ( geneX < 1 || geneY < 0 || geneY >= geneX )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Gene Pair Vector Error"));
      e.setDetails(QObject::tr("Cannot initialize gene vector (%1,%2).").arg(geneX).arg(geneY));
      throw e;
   }
}






Vector::Vector(qint64 index):
   _geneX(1),
   _geneY(0)
{
   // make sure index given is valid
   if ( index < 0 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Gene Pair Vector Error"));
      e.setDetails(QObject::tr("Cannot initialize gene vector from index %1.").arg(index));
      throw e;
   }

   // compute vector from index
   qint64 pos {0};
   while ( pos <= index )
   {
      ++_geneX;
      pos = _geneX * (_geneX - 1) / 2;
   }

   --_geneX;
   pos = _geneX * (_geneX - 1) / 2;

   _geneY = index - pos;
}






qint64 Vector::indent(qint8 cluster) const
{
   // make sure cluster given is valid
   if ( cluster < 0 || cluster >= MAX_CLUSTER_SIZE )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Gene Pair Vector Error"));
      e.setDetails(QObject::tr("Cluster %1 is outside limits.").arg(cluster));
      throw e;
   }

   // compute indent with given cluster and return it
   qint64 index {_geneX * (_geneX - 1) / 2 + _geneY};
   return index * MAX_CLUSTER_SIZE + cluster;
}






void Vector::operator++()
{
   // increment gene y and check it is reaches gene x
   if ( ++_geneY >= _geneX )
   {
      // reset gene y to 0 and increment gene x
      _geneY = 0;
      ++_geneX;
   }
}






Vector Vector::operator++(int)
{
   // save vector value, increment it, and return previous value
   Vector ret {*this};
   ++(*this);
   return ret;
}
