#include "corrpower_workblock.h"



/*!
 * Construct a new block with the given index, starting pairwise index,
 * and pair size.
 *
 * @param index
 * @param start
 * @param size
 */
CorrPowerFilter::WorkBlock::WorkBlock(int index, qint64 start, qint64 size, Pairwise::Index startIndex) :
    EAbstractAnalyticBlock(index),
    _start(start),
    _size(size),
    _startIndex(startIndex)
{
    EDEBUG_FUNC(this,index,startpair,size);
}


/*!
 * Write this block's data to the given data stream.
 *
 * @param stream
 */
void CorrPowerFilter::WorkBlock::write(QDataStream& stream) const
{
    EDEBUG_FUNC(this,&stream);
    stream << _start << _size << _startIndex.getX() << _startIndex.getY();
}



/*!
 * Read this block's data from the given data stream.
 *
 * @param stream
 */
void CorrPowerFilter::WorkBlock::read(QDataStream& stream)
{
    EDEBUG_FUNC(this,&stream);
    qint32 x,y;
    stream >> _start >> _size >> x >> y;
    _startIndex = Pairwise::Index(x,y);
}
