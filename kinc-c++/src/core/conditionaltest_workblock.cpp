#include "conditionaltest_workblock.h"



/*!
 * Implements the interface to create an uninitialized work block at a given index.
 *
 * @param index The given index to create the block at
 */
ConditionalTest::WorkBlock::WorkBlock(int index, qint64 start, qint64 size) :
    EAbstractAnalyticBlock(index),
    _start(start),
    _size(size)
{
    EDEBUG_FUNC(this,index,size);
}



/*!
 * Writes this blocks data to the given data stream.
 *
 * @param stream The data stream that is used to write out data.
 */
void ConditionalTest::WorkBlock::write(QDataStream& stream) const
{
    EDEBUG_FUNC(this,&stream);

    stream << _start << _size;
}



/*!
 * Reads a blocks data from the given data stream.
 *
 * @param stream The data stream that is used to write in data.
 */
void ConditionalTest::WorkBlock::read(QDataStream& stream)
{
    EDEBUG_FUNC(this,&stream);

    stream >> _start >> _size;
}
