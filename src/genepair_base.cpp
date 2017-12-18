#include "genepair_base.h"



using namespace GenePair;






void Base::readData()
{
   // seek to beginning of data object
   if ( !seek(0) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed calling seek() on data object file."));
      throw e;
   }

   // read in all data
   stream() >> _geneSize >> _dataSize >> _pairSize >> _rawPairSize >> _offset;
   readHeader();

   // make sure reading worked
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed reading in data object file."));
      throw e;
   }
}






void Base::newData()
{
   // seek to beginning of data and make sure it worked
   if ( !seek(0) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed calling seek() on data object file."));
      throw e;
   }

   // write out all header information
   stream() << _geneSize << _dataSize << _pairSize << _rawPairSize << _offset;
   writeHeader();

   // make sure writing worked
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed writing to data object file."));
      throw e;
   }
}






const EMetadata& Base::geneNames() const
{
   // get metadata root and make sure genes key exist
   const EMetadata::Map* map {meta().toObject()};
   if ( !map->contains("genes") )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Null Return Reference"));
      e.setDetails(QObject::tr("Requesting reference to gene names when none exists."));
      throw e;
   }

   // return gene names list
   return *(*map)["genes"];
}






void Base::initialize(const EMetadata& geneNames, int dataSize, int offset)
{
   // make sure gene names metadata is an array and is not empty
   if ( !geneNames.isArray() || geneNames.toArray()->isEmpty() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Domain Error"));
      e.setDetails(QObject::tr("Gene names metadata is not an array or empty."));
      throw e;
   }

   // make sure arguments are valid
   if ( dataSize < 1 || offset < 0 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Gene Pair Base Initialization Error"));
      e.setDetails(QObject::tr("An integer argument given is less than the minimum."));
      throw e;
   }

   // get map of metadata root and make copy of gene names
   EMetadata::Map* map {meta().toObject()};
   map->insert("genes",new EMetadata(geneNames));

   // initiailze new data within object
   _geneSize = geneNames.toArray()->size();
   _dataSize = dataSize;
   _offset = offset;
   _pairSize = 0;
   _rawPairSize = 0;
   _lastWrite = -1;
}






void Base::write(Vector index, qint8 cluster)
{
   // make sure this is new data object that can be written to
   if ( _lastWrite == -2 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Gene Pair Base Logical Error"));
      e.setDetails(QObject::tr("Attempting to write data to uninitialized object."));
      throw e;
   }

   // make sure the new gene pair has a higher indent than the previous written so the list of
   // all indents are sorted
   if ( index.indent(cluster) <= _lastWrite )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Gene Pair Base Logical Error"));
      e.setDetails(QObject::tr("Attempting to write indent %1 when last written is %2.")
                   .arg(index.indent(0)).arg(_lastWrite));
      throw e;
   }

   // seek to position for next gene pair and write indent value
   if ( !seek(_headerSize + _offset + _rawPairSize*(_dataSize + _itemHeaderSize)) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed calling seek() on data object file."));
      throw e;
   }
   stream() << index.geneX() << index.geneY() << cluster;

   // make sure writing worked
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed writing to data object file."));
      throw e;
   }

   // increment pair size and set new last index
   ++_rawPairSize;
   _lastWrite = index.indent(cluster);
}






Vector Base::getPair(qint64 index, qint8* cluster) const
{
   // seek to gene pair position and read item header data
   seekPair(index);
   qint32 geneX;
   qint32 geneY;
   stream() >> geneX >> geneY >> *cluster;

   // make sure reading worked
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed reading in gene pair indent."));
      throw e;
   }

   // return gene pair's vector
   return {geneX,geneY};
}






qint64 Base::findPair(qint64 indent, qint64 first, qint64 last) const
{
   // calculate the midway pivot point and seek to it
   qint64 pivot {first + (last - first)/2};
   seekPair(pivot);

   // read in gene pair item header
   qint32 geneX;
   qint32 geneY;
   qint8 cluster;
   stream() >> geneX >> geneY >> cluster;
   Vector index(geneX,geneY);

   // make sure reading worked
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed reading in gene pair indent."));
      throw e;
   }

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






void Base::seekPair(qint64 index) const
{
   // make sure index is within range
   if ( index < 0 || index >= _pairSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Domain Error"));
      e.setDetails(QObject::tr("Attempting to seek to gene pair %1 when total size is %2.")
                   .arg(index).arg(_pairSize));
      throw e;
   }

   // seek to gene pair index requested making sure it worked
   if ( !seek(_headerSize + _offset + index*(_dataSize + _itemHeaderSize)) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed calling seek() on data object file."));
      throw e;
   }
}






void Base::Pair::write(Vector index)
{
   // make sure cluster size of gene pair does not exceed max
   if ( clusterSize() > Vector::_maxClusterSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Gene Pair Logical Error"));
      e.setDetails(QObject::tr("Cannot write gene pair with cluster size of %1 exceeding the max of"
                               " %2.").arg(clusterSize()).arg(Vector::_maxClusterSize));
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






void Base::Pair::read(Vector index) const
{
   // clear any existing clusters
   clearClusters();

   // attempt to find vector within data object
   qint64 index_;
   if ( _cMatrix->_rawPairSize > 0
        && (index_ = _cMatrix->findPair(index.indent(0),0,_cMatrix->_rawPairSize - 1)) != -1 )
   {
      // gene pair with vector found, read in all clusters
      _nextIndex = index_;
      readNext();
   }
}






void Base::Pair::readFirst() const
{
   // reset next index to first gene pair and read it
   _nextIndex = 0;
   readNext();
}






void Base::Pair::readNext() const
{
   // make sure read next index is not already at end of data object
   if ( _nextIndex < _cMatrix->_rawPairSize )
   {
      // clear any existing clusters
      clearClusters();

      // get to first cluster
      qint8 cluster;
      Vector vector {_cMatrix->getPair(_nextIndex++,&cluster)};

      // make sure this is cluster 0
      if ( cluster != 0 )
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(QObject::tr("File IO Error"));
         e.setDetails(QObject::tr("Reading gene pair failed because first cluster is not 0."));
         throw e;
      }

      // add first cluster, read it in, and save what vector this gene pair is
      addCluster();
      readCluster(_cMatrix->stream(),0);
      _vector = vector;

      // read in remaining clusters for gene pair
      qint8 count {1};
      while ( _nextIndex < _cMatrix->_rawPairSize )
      {
         // get next gene pair cluster
         _cMatrix->getPair(_nextIndex++,&cluster);

         // if cluster is zero this is the next gene pair so break from loop
         if ( cluster == 0 )
         {
            --_nextIndex;
            break;
         }

         // make sure max cluster size has not been exceeded
         if ( ++count > Vector::_maxClusterSize )
         {
            E_MAKE_EXCEPTION(e);
            e.setTitle(QObject::tr("File IO Error"));
            e.setDetails(QObject::tr("Reading gene pair failed because it exceeds the max number of"
                                     " clusters."));
            throw e;
         }

         // add new cluster and read it in
         addCluster();
         readCluster(_cMatrix->stream(),cluster);
      }
   }
}
