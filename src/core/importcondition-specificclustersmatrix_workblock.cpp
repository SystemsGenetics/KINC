#include "importcondition-specificclustersmatrix_workblock.h"
//





/*!
*  Implements the interface to create an uninitialized work block at a given index.
*
*  @param index The given index to create the block at
*/
importCSCM::WorkBlock::WorkBlock(int index, Pairwise::Index start, qint64 startpair, qint64 size) :
    EAbstractAnalyticBlock(index),
    _start(start),
    _startpair(startpair),
    _size(size)
{
    EDEBUG_FUNC(this,index,startpair,size);
}





/*!
*  Writes this blocks data to the given data stream.
*
*  @param stream The data stream that is used to write out data.
*/
void importCSCM::WorkBlock::write(QDataStream& stream) const
{
    EDEBUG_FUNC(this,&stream);

    stream << _startpair << _size;
}





/*!
*  Reads a blocks data from the given data stream.
*
*  @param stream The data stream that is used to write in data.
*/
void importCSCM::WorkBlock::read(QDataStream& stream)
{
    EDEBUG_FUNC(this,&stream);

    stream >> _startpair >> _size;
}
