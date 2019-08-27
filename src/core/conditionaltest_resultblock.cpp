#include "conditionaltest_resultblock.h"
//




/*!
*  Implements an interface to create a result block object.
*/
ConditionalTest::ResultBlock::ResultBlock(int index) : EAbstractAnalyticBlock(index)
{
    EDEBUG_FUNC(this,index);
}





/*!
*  Implements an interface to create a result block object.
*/
ConditionalTest::ResultBlock::ResultBlock(int index, int numTests, qint64 start) :
    EAbstractAnalyticBlock(index),
    _numTests(numTests),
    _start(start)
{
    EDEBUG_FUNC(this,index,numTests,start);
}






/*!
 * Append a pair to the result block's list of pairs.
 *
 * @param pair
 */
void ConditionalTest::ResultBlock::append(const CSCMPair& pair)
{
   EDEBUG_FUNC(this,&pair);

   _pairs.append(pair);
}





/*!
*  Implements an interface to write the result block header into the file stream.
*/
void ConditionalTest::ResultBlock::write(QDataStream& stream) const
{
    EDEBUG_FUNC(this,&stream);
    for(auto cell : _pairs)
    {
        stream << cell.pValues.size();
        for(int i = 0; i < cell.pValues.size(); i++)
        {
            stream << cell.pValues.at(i);
        }
    }
}





/*!
*  Implements an interface to read the result block header from the file stream.
*/
void ConditionalTest::ResultBlock::read(QDataStream& stream)
{
    EDEBUG_FUNC(this,&stream);
    int size = 0;
    for(int i = 0; i < _pairs.size(); i++)
    {
        stream >> size;
        for(int j = 0; j < size; j++)
        {
            stream << _pairs[i].pValues[j];
        }
    }
}
