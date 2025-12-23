#include "similarity_resultblock.h"



/*!
 * Construct a new block with the given index.
 *
 * @param index
 */
Similarity::ResultBlock::ResultBlock(int index):
    EAbstractAnalyticBlock(index)
{
    EDEBUG_FUNC(this,index);
}



/*!
 * Append a pair to the result block's list of pairs.
 *
 * @param pair
 */
void Similarity::ResultBlock::append(const Pair& pair)
{
    EDEBUG_FUNC(this,&pair);

    _pairs.append(pair);
}



/*!
 * Write this block's data to the given data stream.
 *
 * @param stream
 */
void Similarity::ResultBlock::write(QDataStream& stream) const
{
    EDEBUG_FUNC(this,&stream);

    stream << _pairs.size();

    for ( auto& pair : _pairs )
    {
        stream << pair.index.getX();
        stream << pair.index.getY();
        stream << pair.labels;
        stream << pair.correlations;
    }
}



/*!
 * Read this block's data from the given data stream.
 *
 * @param stream
 */
void Similarity::ResultBlock::read(QDataStream& stream)
{
    EDEBUG_FUNC(this,&stream);

    int size;
    stream >> size;

    _pairs.resize(size);

    for ( auto& pair : _pairs )
    {
        qint32 x, y;
        stream >> x;
        stream >> y;

        pair.index = Pairwise::Index(x, y);

        stream >> pair.labels;
        stream >> pair.correlations;
    }
}
