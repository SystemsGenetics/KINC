#include "genepair_vector.h"




using namespace GenePair;






Vector::Vector(qint64 geneX, qint64 geneY, qint64 cluster):
   _geneX(geneX),
   _geneY(geneY),
   _cluster(cluster)
{
   if ( geneX < 1 || geneY < 0 || cluster < 0 || geneY >= geneX || cluster >= _maxClusterSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Gene Pair Vector Error"));
      e.setDetails(QObject::tr("Cannot initialize gene vector (%1,%2,%3).").arg(geneX).arg(geneY)
                   .arg(cluster));
      throw e;
   }
}






void Vector::indent(qint64 indent, qint64 geneX)
{
   qint64 cluster = indent%_maxClusterSize;
   qint64 geneY = indent/_maxClusterSize - geneX*(geneX - 1);
   if ( geneX < 1 || geneY < 0 || cluster < 0 || geneY >= geneX || cluster >= _maxClusterSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Gene Pair Vector Error"));
      e.setDetails(QObject::tr("Cannot initialize gene vector (%1,%2,%3).").arg(geneX).arg(geneY)
                   .arg(cluster));
      throw e;
   }
   _geneX = geneX;
   _geneY = geneY;
   _cluster = cluster;
}






void Vector::operator++()
{
   if ( ++_cluster >= _maxClusterSize )
   {
      _cluster = 0;
      if ( ++_geneY >= _geneX )
      {
         _geneY = 0;
         ++_geneX;
      }
   }
}






Vector Vector::operator++(int)
{
   Vector ret {*this};
   ++(*this);
   return ret;
}
