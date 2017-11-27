#include "correlationbase.h"






void CorrelationBase::initialize(int geneSize, int dataSize, int offset)
{
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
   if ( y >= x )
   {
      ;//ERROR
   }
}






void CorrelationBase::Iterator::operator--()
{
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
   if ( ++_y >= _x )
   {
      _y = 0;
      ++_x;
   }
}
