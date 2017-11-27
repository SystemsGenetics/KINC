#include "correlationbase.h"






void CorrelationBase::finish()
{
   // seek to beginning of data object
   if ( !EAbstractData::seek(0) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(QObject::tr("File IO Error"));
      e.setDetails(QObject::tr("Failed calling seek() on data object file."));
      throw e;
   }

   // write out all data
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






void CorrelationBase::write(Iterator correlation)
{
}






void CorrelationBase::read()
{
}






bool CorrelationBase::seek(CorrelationBase::Iterator correlation) const
{
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
