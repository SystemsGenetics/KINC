#include "conditionaltest_workblock.h"



/*!
 * Implements the interface to create an uninitialized work block at a given index.
 *
 * @param index The given index to create the block at
 */
ConditionalTest::WorkBlock::WorkBlock(int index, Pairwise::Index start, qint64 startpair, qint64 size) :
    EAbstractAnalyticBlock(index),
    _start(start),
    _startpair(startpair),
    _size(size)
{
    EDEBUG_FUNC(this,index,startpair,size);
}



/*!
 * Writes this blocks data to the given data stream.
 *
 * @param stream The data stream that is used to write out data.
 */
void ConditionalTest::WorkBlock::write(QDataStream& stream) const
{
    EDEBUG_FUNC(this,&stream);

    stream << _startpair << _size;
}



/*!
 * Reads a blocks data from the given data stream.
 *
 * @param stream The data stream that is used to write in data.
 */
void ConditionalTest::WorkBlock::read(QDataStream& stream)
{
    EDEBUG_FUNC(this,&stream);

    stream >> _startpair >> _size;
}
