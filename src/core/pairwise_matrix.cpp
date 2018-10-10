#include "pairwise_matrix.h"



using namespace Pairwise;






/*!
 * Return the index of the first byte in this data object after the end of
 * the data section. Defined as the size of the header and sub-header plus the
 * total size of all pairs.
 */
qint64 Matrix::dataEnd() const
{
   EDEBUG_FUNC(this);

   return _headerSize + _subHeaderSize + _clusterSize * (_dataSize + _itemHeaderSize);
}






/*!
 * Read in the data of an existing data object that was just opened.
 */
void Matrix::readData()
{
   EDEBUG_FUNC(this);

   // seek to the beginning of the data
   seek(0);

   // read the header
   stream() >> _geneSize >> _maxClusterSize >> _dataSize >> _pairSize >> _clusterSize >> _subHeaderSize;

   // read the sub-header
   readHeader();
}






/*!
 * Initialize this data object's data to a null state.
 */
void Matrix::writeNewData()
{
   EDEBUG_FUNC(this);

   // initialize metadata
   setMeta(EMetadata(EMetadata::Object));

   // seek to the beginning of the data
   seek(0);

   // write the header
   stream() << _geneSize << _maxClusterSize << _dataSize << _pairSize << _clusterSize << _subHeaderSize;

   // write the sub-header
   writeHeader();
}






/*!
 * Finalize this data object's data after the analytic that created it has
 * finished giving it new data.
 */
void Matrix::finish()
{
   EDEBUG_FUNC(this);

   // seek to the beginning of the data
   seek(0);

   // write the header
   stream() << _geneSize << _maxClusterSize << _dataSize << _pairSize << _clusterSize << _subHeaderSize;

   // write the sub-header
   writeHeader();
}






/*!
 * Return the list of gene names in this pairwise matrix.
 */
EMetadata Matrix::geneNames() const
{
   EDEBUG_FUNC(this);

   return meta().toObject().at("genes");
}






/*!
 * Initialize this pairwise matrix with a list of gene names, the max cluster
 * size, the pairwise data size, and the sub-header size.
 *
 * @param geneNames
 * @param maxClusterSize
 * @param dataSize
 * @param subHeaderSize
 */
void Matrix::initialize(const EMetadata& geneNames, int maxClusterSize, int dataSize, int subHeaderSize)
{
   EDEBUG_FUNC(this,geneNames,maxClusterSize,dataSize,subHeaderSize);

   // make sure gene names metadata is an array and is not empty
   if ( !geneNames.isArray() || geneNames.toArray().isEmpty() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Gene names metadata is not an array or is empty."));
      throw e;
   }

   // make sure arguments are valid
   if ( maxClusterSize < 1 || dataSize < 1 || subHeaderSize < 0 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Pairwise Matrix Initialization Error"));
      e.setDetails(tr("An integer argument given is less than the minimum."));
      throw e;
   }

   if ( maxClusterSize > Index::MAX_CLUSTER_SIZE )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Pairwise Matrix Initialization Error"));
      e.setDetails(tr("The max cluster size cannot be larger than MAX_CLUSTER_SIZE."));
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
   _subHeaderSize = subHeaderSize;
   _pairSize = 0;
   _clusterSize = 0;
   _lastWrite = -1;
}






/*!
 * Write the header of a new pair given a pairwise index and cluster index.
 *
 * @param index
 * @param cluster
 */
void Matrix::write(const Index& index, qint8 cluster)
{
   EDEBUG_FUNC(this,index,cluster);

   // make sure this is new data object that can be written to
   if ( _lastWrite == -2 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Pairwise Matrix Logical Error"));
      e.setDetails(tr("Attempting to write data to uninitialized object."));
      throw e;
   }

   // make sure the new pair has a higher indent than the previous written so the list of
   // all indents are sorted
   if ( index.indent(cluster) <= _lastWrite )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Pairwise Matrix Logical Error"));
      e.setDetails(tr("Attempting to write indent %1 when last written is %2.")
                   .arg(index.indent(0)).arg(_lastWrite));
      throw e;
   }

   // seek to position for next pair and write indent value
   seek(_headerSize + _subHeaderSize + _clusterSize * (_dataSize + _itemHeaderSize));
   stream() << index.getX() << index.getY() << cluster;

   // increment cluster size and set new last index
   ++_clusterSize;
   _lastWrite = index.indent(cluster);
}






/*!
 * Get a pair at the given index in the data object file and return the
 * pairwise index and cluster index of that pair.
 *
 * @param index
 * @param cluster
 */
Index Matrix::getPair(qint64 index, qint8* cluster) const
{
   EDEBUG_FUNC(this,index,cluster);

   // seek to index and read item header data
   seekPair(index);
   qint32 geneX;
   qint32 geneY;
   stream() >> geneX >> geneY >> *cluster;

   // return pairwise index
   return {geneX,geneY};
}






/*!
 * Find a pair with a given indent value using binary search.
 *
 * @param indent
 * @param first
 * @param last
 */
qint64 Matrix::findPair(qint64 indent, qint64 first, qint64 last) const
{
   EDEBUG_FUNC(this,indent,first,last);

   // calculate the midway pivot point and seek to it
   qint64 pivot {first + (last - first)/2};
   seekPair(pivot);

   // read in pairwise item header
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

   // else no pair with given indent exists so return failure
   else
   {
      return -1;
   }
}






/*!
 * Seek to the pair at the given index in the data object file.
 *
 * @param index
 */
void Matrix::seekPair(qint64 index) const
{
   EDEBUG_FUNC(this,index);

   // make sure index is within range
   if ( index < 0 || index >= _clusterSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to seek to cluster index %1 when total size is %2.")
                   .arg(index).arg(_clusterSize));
      throw e;
   }

   // seek to the specified index
   seek(_headerSize + _subHeaderSize + index * (_dataSize + _itemHeaderSize));
}
