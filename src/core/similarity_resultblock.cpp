#include "similarity_resultblock.h"



/*!
 * Construct a new block with the given index and starting pairwise index.
 *
 * @param index
 * @param start
 */
Similarity::ResultBlock::ResultBlock(int index, qint64 start):
   EAbstractAnalyticBlock(index),
   _start(start)
{
   EDEBUG_FUNC(this,index,start);
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

   stream << _start;
   stream << _pairs.size();

   for ( auto& pair : _pairs )
   {
      stream << pair.K;
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

   stream >> _start;

   int size;
   stream >> size;

   _pairs.resize(size);

   for ( auto& pair : _pairs )
   {
      stream >> pair.K;
      stream >> pair.labels;
      stream >> pair.correlations;
   }
}
