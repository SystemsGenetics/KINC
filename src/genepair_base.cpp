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
   stream() >> _geneSize >> _dataSize >> _pairSize >> _offset;

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
   stream() << _geneSize << _dataSize << _pairSize << _offset;

   // make sure writing worked
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed writing to data object file."));
      throw e;
   }
}






void Base::prepare(bool preAllocate)
{
   // check to see if pre-allocation is requested
   if ( preAllocate )
   {
      // seek to beginning of data and make sure it worked
      if ( !seek(0) )
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(QObject::tr("File IO Error"));
         e.setDetails(QObject::tr("Failed calling seek() on data object file."));
         throw e;
      }

      // allocate total size needed in data and make sure it worked
      if ( !allocate(getDataEnd()) )
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(QObject::tr("File IO Error"));
         e.setDetails(QObject::tr("Failed allocating space in data object file."));
         throw e;
      }
   }
}






void Base::initialize(int geneSize, int dataSize, int offset)
{
   // make sure arguments are valid
   if ( geneSize < 1 || dataSize < 1 || offset < 0 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Correlation Base Initialization Error"));
      e.setDetails(QObject::tr("An integer argument given is less than the minimum."));
      throw e;
   }

   // initiailze new data within object
   _geneSize = geneSize;
   _dataSize = dataSize;
   _pairSize = 0;
   _offset = offset;
   _lastWrite = Vector();
}






void Base::write(Vector index)
{
   // make sure the new correlation has a higher indent than the previous written so the list of
   // all indents are sorted
   if ( index.indent() >= _lastWrite.indent() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Correlation Base Logical Error"));
      e.setDetails(QObject::tr("Attempting to write indent %1 when last written is %2.")
                   .arg(index.indent()).arg(_lastWrite.indent()));
      throw e;
   }

   // seek to position for next correlation and write indent value
   seekPair(_pairSize);
   stream() << index.indent();

   // make sure writing worked
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed writing to data object file."));
      throw e;
   }

   // increment pair size and set new last index
   ++_pairSize;
   _lastWrite = index;
}






qint64 Base::getPair(qint64 index) const
{
   // seek to correlation position and read indent data
   seekPair(index);
   qint64 indent;
   stream() >> indent;

   // make sure reading worked and return indent
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed reading in gene pair indent."));
      throw e;
   }
   return indent;
}






qint64 Base::findPair(qint64 indent, qint64 first, qint64 last) const
{
   // calculate the midway pivot point and seek to it
   qint64 pivot {first + (last - first)/2};
   seekPair(pivot);

   // read in indent value of correlation
   qint64 readIndent;
   stream() >> readIndent;

   // make sure reading worked
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed reading in gene pair indent."));
      throw e;
   }

   // if indent values match return index
   if ( readIndent == indent )
   {
      return pivot;
   }

   // else if first and last are not the same divide and conquer
   else if ( first != last )
   {
      // if indent is less than pivot divide lower half
      if ( readIndent > indent )
      {
         // if pivot is first add one so pivot is not less than first when passed
         if ( pivot == first )
         {
            ++pivot;
         }
         return findPair(indent,first,pivot - 1);
      }

      // else it is greater so divide upper half
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

   // else no correlation with given indent exists so return failure
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

   // seek to correlation index requested making sure it worked
   if ( !seek(_headerSize + _offset + index*(_dataSize + sizeof(qint64))) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed calling seek() on data object file."));
      throw e;
   }
}
