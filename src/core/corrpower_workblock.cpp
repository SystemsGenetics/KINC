#include "corrpower_workblock.h"



/*!
 * Construct a new block with the given index, starting pairwise index,
 * and pair size.
 *
 * @param index
 * @param start
 * @param size
 */
CorrPowerFilter::WorkBlock::WorkBlock(int index, qint64 start, qint64 size):
   EAbstractAnalyticBlock(index),
   _start(start),
   _size(size)
{
   EDEBUG_FUNC(this,index,start,size);
}



/*!
 * Write this block's data to the given data stream.
 *
 * @param stream
 */
void CorrPowerFilter::WorkBlock::write(QDataStream& stream) const
{
   EDEBUG_FUNC(this,&stream);

   stream << _start << _size;
}



/*!
 * Read this block's data from the given data stream.
 *
 * @param stream
 */
void CorrPowerFilter::WorkBlock::read(QDataStream& stream)
{
   EDEBUG_FUNC(this,&stream);

   stream >> _start >> _size;
}
