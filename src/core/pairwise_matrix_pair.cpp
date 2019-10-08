#include "pairwise_matrix_pair.h"



using namespace Pairwise;






/*!
 * Write the iterator's pairwise data to the data object file with the given
 * pairwise index.
 *
 * @param index
 */
void Matrix::Pair::write(const Index& index)
{
   EDEBUG_FUNC(this,&index);

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
   for ( qint8 i = 0; i < clusterSize(); ++i )
   {
      _matrix->write(index,i);
      writeCluster(_matrix->stream(),i);
   }

   // increment pair size of data object
   ++(_matrix->_pairSize);
}






/*!
 * Read the pair with the given pairwise index from the data object file.
 *
 * @param index
 */
void Matrix::Pair::read(const Index& index) const
{
   EDEBUG_FUNC(this,&index);

   // clear any existing clusters
   clearClusters();

   // attempt to find cluster index within data object
   qint64 clusterIndex {_cMatrix->findPair(index.indent(0), 0, _cMatrix->_clusterSize - 1)};

   if ( _cMatrix->_clusterSize > 0 && clusterIndex != -1 )
   {
      // pair found, read in all clusters
      _rawIndex = clusterIndex;
      readNext();
   }
}






/*!
 * Read the next pair in the data object file.
 */
void Matrix::Pair::readNext() const
{
   EDEBUG_FUNC(this);

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
