#include "genepair_vector.h"




using namespace GenePair;






Vector::Vector(qint32 geneX, qint32 geneY):
   _geneX(geneX),
   _geneY(geneY)
{
   if ( geneX < 1 || geneY < 0 || geneY >= geneX )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Gene Pair Vector Error"));
      e.setDetails(QObject::tr("Cannot initialize gene vector (%1,%2).").arg(geneX).arg(geneY));
      throw e;
   }
}






qint64 Vector::indent(qint8 cluster) const
{
   if ( cluster < 0 || cluster >= _maxClusterSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Gene Pair Vector Error"));
      e.setDetails(QObject::tr("Cluster %1 is outside limits.").arg(cluster));
      throw e;
   }
   qint64 ret {_geneX};
   ret = ret*(ret - 1)/2;
   ret += _geneY;
   ret = ret*_maxClusterSize + cluster;
   return ret;
}






void Vector::operator++()
{
   if ( ++_geneY >= _geneX )
   {
      _geneY = 0;
      ++_geneX;
   }
}






Vector Vector::operator++(int)
{
   Vector ret {*this};
   ++(*this);
   return ret;
}
