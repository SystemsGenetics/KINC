#include "conditionaltest_resultblock.h"



/*!
 * Create a result block object.
 */
ConditionalTest::ResultBlock::ResultBlock(int index) : EAbstractAnalyticBlock(index)
{
    EDEBUG_FUNC(this,index);
}



/*!
 * Append a pair to the result block's list of pairs.
 *
 * @param pair
 */
void ConditionalTest::ResultBlock::append(const Pair& pair)
{
    EDEBUG_FUNC(this,&pair);

    _pairs.append(pair);
}



/*!
 * Write the result block header into the file stream.
 */
void ConditionalTest::ResultBlock::write(QDataStream& stream) const
{
    EDEBUG_FUNC(this,&stream);

    stream << _pairs.size();

    for ( auto& pair : _pairs )
    {
        // write pairwise index
        stream << pair.index.getX();
        stream << pair.index.getY();

        // write p values and r^2 values
        stream << pair.pValues.size();

        for ( int k = 0; k < pair.pValues.size(); k++ )
        {
            stream << pair.pValues.at(k);
            stream << pair.r2.at(k);
        }
    }
}



/*!
 * Read the result block header from the file stream.
 */
void ConditionalTest::ResultBlock::read(QDataStream& stream)
{
    EDEBUG_FUNC(this,&stream);

    int size;
    stream >> size;

    _pairs.resize(size);

    for ( auto& pair : _pairs )
    {
        // read pairwise index
        qint32 x, y;
        stream >> x;
        stream >> y;

        pair.index = Pairwise::Index(x, y);

        // read p values and r^2 values
        stream >> size;

        pair.pValues.resize(size);
        pair.r2.resize(size);

        for ( int k = 0; k < pair.pValues.size(); k++ )
        {
            stream >> pair.pValues[k];
            stream >> pair.r2[k];
        }
    }
}
