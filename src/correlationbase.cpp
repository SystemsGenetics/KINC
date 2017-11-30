#include "correlationbase.h"






void CorrelationBase::readData()
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
   stream() >> _geneSize >> _dataSize >> _correlationSize >> _offset;

   // make sure reading worked
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed reading in data object file."));
      throw e;
   }
}






void CorrelationBase::newData()
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
   stream() << _geneSize << _dataSize << _correlationSize << _offset;

   // make sure writing worked
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed writing to data object file."));
      throw e;
   }
}






void CorrelationBase::prepare(bool preAllocate)
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






void CorrelationBase::initialize(int geneSize, int dataSize, int offset)
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
   _correlationSize = 0;
   _offset = offset;
   _lastWrite = Iterator();

   // seek to end of header data for base class
   if ( !EAbstractData::seek(_headerSize) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed calling seek() on data object file."));
      throw e;
   }
}






void CorrelationBase::write(Iterator index)
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
   seekCorrelation(_correlationSize++);
   stream() << index.indent();

   // make sure writing worked and set new last written correlation
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed writing to data object file."));
      throw e;
   }
   _lastWrite = index;
}






bool CorrelationBase::findCorrelation(qint64 indent, int first, int last) const
{
   // calculate the midway pivot point and seek to it
   int pivot {first + (last - first)/2};
   seekCorrelation(pivot);

   // read in indent value of correlation
   qint64 readIndent;
   stream() >> readIndent;

   // make sure reading worked
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed reading in correlation indent."));
      throw e;
   }

   // if indent values match return true
   if ( readIndent == indent )
   {
      return true;
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
         return findCorrelation(indent,first,pivot - 1);
      }

      // else it is greater so divide upper half
      else
      {
         // if pivot is last decrement one so pivot is not greater than last when passed
         if ( pivot == last )
         {
            --pivot;
         }
         return findCorrelation(indent,pivot + 1,last);
      }
   }

   // else no correlation with given indent exists so return false
   else
   {
      return false;
   }
}






void CorrelationBase::seekCorrelation(int index) const
{
   // seek to correlation index requested making sure it worked
   if ( !seek(_headerSize + _offset + index*(_dataSize + sizeof(qint64))) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed calling seek() on data object file."));
      throw e;
   }
}






CorrelationBase::Iterator::Iterator(int x, int y):
   _x(x),
   _y(y)
{
   // make sure the gene pair is valid
   if ( x < 1 || y < 0 || y >= x )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("Correlation Iterator Error"));
      e.setDetails(QObject::tr("Cannot initialize correlation pair (%1,%2).").arg(x).arg(y));
      throw e;
   }
}






void CorrelationBase::Iterator::operator--()
{
   // decrement the gene pair
   if ( --_y < 0 )
   {
      if ( --_x < 1 )
      {
         _x = 1;
      }
      _y = _x - 1;
   }
}






void CorrelationBase::Iterator::operator++()
{
   // increment the gene pair
   if ( ++_y >= _x )
   {
      _y = 0;
      ++_x;
   }
}
