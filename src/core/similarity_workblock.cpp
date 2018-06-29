#include "similarity_workblock.h"






Similarity::WorkBlock::WorkBlock(int index, qint64 start, qint64 size):
   EAbstractAnalytic::Block(index),
   _start(start),
   _size(size)
{
}






void Similarity::WorkBlock::write(QDataStream& stream) const
{
   stream << _start << _size;
}






void Similarity::WorkBlock::read(QDataStream& stream)
{
   stream >> _start >> _size;
}
