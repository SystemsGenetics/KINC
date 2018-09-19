#include "similarity_workblock.h"






/*!
 * Construct a new block with the given index, starting pairwise index,
 * and pair size.
 *
 * @param index
 * @param start
 * @param size
 */
Similarity::WorkBlock::WorkBlock(int index, qint64 start, qint64 size):
   EAbstractAnalytic::Block(index),
   _start(start),
   _size(size)
{
}






/*!
 * Write this block's data to the given data stream.
 *
 * @param stream
 */
void Similarity::WorkBlock::write(QDataStream& stream) const
{
   stream << _start << _size;
}






/*!
 * Read this block's data from the given data stream.
 *
 * @param stream
 */
void Similarity::WorkBlock::read(QDataStream& stream)
{
   stream >> _start >> _size;
}
