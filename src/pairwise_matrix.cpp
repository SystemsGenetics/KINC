#include "pairwise_matrix.h"



using namespace Pairwise;






qint64 Matrix::dataEnd() const
{
   return _headerSize + _offset + _clusterSize * (_dataSize + _itemHeaderSize);
}






void Matrix::readData()
{
   // read header
   seek(0);
   stream() >> _geneSize >> _maxClusterSize >> _dataSize >> _pairSize >> _clusterSize >> _offset;
   readHeader();
}






void Matrix::writeNewData()
{
   // initialize metadata
   setMeta(EMetadata(EMetadata::Object));

   // initialize header
   seek(0);
   stream() << _geneSize << _maxClusterSize << _dataSize << _pairSize << _clusterSize << _offset;
   writeHeader();
}






void Matrix::finish()
{
   // initialize header
   seek(0);
   stream() << _geneSize << _maxClusterSize << _dataSize << _pairSize << _clusterSize << _offset;
   writeHeader();
}






EMetadata Matrix::geneNames() const
{
   return meta().toObject().at("genes");
}






void Matrix::initialize(const EMetadata& geneNames, int maxClusterSize, int dataSize, int offset)
{
   // make sure gene names metadata is an array and is not empty
   if ( !geneNames.isArray() || geneNames.toArray().isEmpty() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Gene names metadata is not an array or is empty."));
      throw e;
   }

   // make sure arguments are valid
   if ( maxClusterSize < 1 || dataSize < 1 || offset < 0 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Pairwise Matrix Initialization Error"));
      e.setDetails(tr("An integer argument given is less than the minimum."));
      throw e;
   }

   // save gene names to metadata
   EMetaObject metaObject {meta().toObject()};
   metaObject.insert("genes", geneNames);
   setMeta(metaObject);

   // initiailze new data within object
   _geneSize = geneNames.toArray().size();
   _maxClusterSize = maxClusterSize;
   _dataSize = dataSize;
   _offset = offset;
   _pairSize = 0;
   _clusterSize = 0;
   _lastWrite = -1;
}






void Matrix::write(Index index, qint8 cluster)
{
   // make sure this is new data object that can be written to
   if ( _lastWrite == -2 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Pairwise Matrix Logical Error"));
      e.setDetails(tr("Attempting to write data to uninitialized object."));
      throw e;
   }

   // make sure the new gene pair has a higher indent than the previous written so the list of
   // all indents are sorted
   if ( index.indent(cluster) <= _lastWrite )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Pairwise Matrix Logical Error"));
      e.setDetails(tr("Attempting to write indent %1 when last written is %2.")
                   .arg(index.indent(0)).arg(_lastWrite));
      throw e;
   }

   // seek to position for next gene pair and write indent value
   seek(_headerSize + _offset + _clusterSize * (_dataSize + _itemHeaderSize));
   stream() << index.getX() << index.getY() << cluster;

   // increment cluster size and set new last index
   ++_clusterSize;
   _lastWrite = index.indent(cluster);
}






Index Matrix::getPair(qint64 index, qint8* cluster) const
{
   // seek to gene pair position and read item header data
   seekPair(index);
   qint32 geneX;
   qint32 geneY;
   stream() >> geneX >> geneY >> *cluster;

   // return gene pair index
   return {geneX,geneY};
}






qint64 Matrix::findPair(qint64 indent, qint64 first, qint64 last) const
{
   // calculate the midway pivot point and seek to it
   qint64 pivot {first + (last - first)/2};
   seekPair(pivot);

   // read in gene pair item header
   qint32 geneX;
   qint32 geneY;
   qint8 cluster;
   stream() >> geneX >> geneY >> cluster;
   Index index(geneX,geneY);

   // if indent values match return index
   if ( index.indent(cluster) == indent )
   {
      return pivot;
   }

   // else if the list of indexes is greater than 1 divide it in half and figure out which side
   // to divide and conquer depending on the value of the pivot
   else if ( first != last )
   {
      if ( index.indent(cluster) > indent )
      {
         // if pivot is first add one so pivot is not less than first when passed
         if ( pivot == first )
         {
            ++pivot;
         }
         return findPair(indent,first,pivot - 1);
      }
      else
      {
         // if pivot is last decrement one so pivot is not greater than last when passed
         if ( pivot == last )
         {
            --pivot;
         }
         return findPair(indent,pivot + 1,last);
      }
   }

   // else no gene pair with given indent exists so return failure
   else
   {
      return -1;
   }
}






void Matrix::seekPair(qint64 index) const
{
   // make sure index is within range
   if ( index < 0 || index >= _clusterSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to seek to cluster index %1 when total size is %2.")
                   .arg(index).arg(_clusterSize));
      throw e;
   }

   // seek to gene pair index requested making sure it worked
   seek(_headerSize + _offset + index * (_dataSize + _itemHeaderSize));
}






void Matrix::Pair::write(Index index)
{
   // make sure cluster size of gene pair does not exceed max
   if ( clusterSize() > Index::MAX_CLUSTER_SIZE )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Pairwise Logical Error"));
      e.setDetails(tr("Cannot write gene pair with cluster size of %1 exceeding the max of"
                               " %2.").arg(clusterSize()).arg(Index::MAX_CLUSTER_SIZE));
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
      // gene pair found, read in all clusters
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
         e.setDetails(tr("Reading gene pair failed because first cluster is not 0."));
         throw e;
      }

      // add first cluster, read it in, and save pairwise index
      addCluster();
      readCluster(_cMatrix->stream(),0);
      _index = index;

      // read in remaining clusters for gene pair
      qint8 count {1};
      while ( _rawIndex < _cMatrix->_clusterSize )
      {
         // get next gene pair cluster
         _cMatrix->getPair(_rawIndex++,&cluster);

         // if cluster is zero this is the next gene pair so break from loop
         if ( cluster == 0 )
         {
            --_rawIndex;
            break;
         }

         // make sure max cluster size has not been exceeded
         if ( ++count > Index::MAX_CLUSTER_SIZE )
         {
            E_MAKE_EXCEPTION(e);
            e.setTitle(tr("File IO Error"));
            e.setDetails(tr("Reading gene pair failed because it exceeds the max number of"
                                     " clusters."));
            throw e;
         }

         // add new cluster and read it in
         addCluster();
         readCluster(_cMatrix->stream(),cluster);
      }
   }
}
