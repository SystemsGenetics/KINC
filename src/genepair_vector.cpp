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
