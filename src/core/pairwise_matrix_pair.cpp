#include "pairwise_matrix_pair.h"



using namespace Pairwise;






void Matrix::Pair::write(Index index)
{
   // make sure cluster size of pair does not exceed max
   if ( clusterSize() > _matrix->_maxClusterSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Pairwise Logical Error"));
      e.setDetails(tr("Cannot write pair with cluster size %1 exceeding the max of %2.")
         .arg(clusterSize())
         .arg(_matrix->_maxClusterSize));
      throw e;
   }

   // go through each cluster and write it to data object
   for (int i = 0; i < clusterSize() ;++i)
   {
      _matrix->write(index,i);
      writeCluster(_matrix->stream(),i);
   }

   // increment pair size of data object
   ++(_matrix->_pairSize);
}






void Matrix::Pair::read(Index index) const
{
   // clear any existing clusters
   clearClusters();

   // attempt to find cluster index within data object
   qint64 clusterIndex;
   if ( _cMatrix->_clusterSize > 0
        && (clusterIndex = _cMatrix->findPair(index.indent(0),0,_cMatrix->_clusterSize - 1)) != -1 )
   {
      // pair found, read in all clusters
      _rawIndex = clusterIndex;
      readNext();
   }
}






void Matrix::Pair::readNext() const
{
   // make sure read next index is not already at end of data object
   if ( _rawIndex < _cMatrix->_clusterSize )
   {
      // clear any existing clusters
      clearClusters();

      // get to first cluster
      qint8 cluster;
      Index index {_cMatrix->getPair(_rawIndex++,&cluster)};

      // make sure this is cluster 0
      if ( cluster != 0 )
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("File IO Error"));
         e.setDetails(tr("Reading pair failed because first cluster is not 0."));
         throw e;
      }

      // add first cluster, read it in, and save pairwise index
      addCluster();
      readCluster(_cMatrix->stream(),0);
      _index = index;

      // read in remaining clusters for pair
      qint8 count {1};
      while ( _rawIndex < _cMatrix->_clusterSize )
      {
         // get next pair cluster
         _cMatrix->getPair(_rawIndex++,&cluster);

         // if cluster is zero this is the next pair so break from loop
         if ( cluster == 0 )
         {
            --_rawIndex;
            break;
         }

         // make sure max cluster size has not been exceeded
         if ( ++count > _cMatrix->_maxClusterSize )
         {
            E_MAKE_EXCEPTION(e);
            e.setTitle(tr("Pairwise Logical Error"));
            e.setDetails(tr("Cannot read pair with cluster size %1 exceeding the max of %2.")
               .arg(count)
               .arg(_matrix->_maxClusterSize));
            throw e;
         }

         // add new cluster and read it in
         addCluster();
         readCluster(_cMatrix->stream(),cluster);
      }
   }
}