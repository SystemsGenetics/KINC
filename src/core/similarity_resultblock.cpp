#include "similarity_resultblock.h"






/*!
 * Construct a new block with the given index and starting pairwise index.
 *
 * @param index
 * @param start
 */
Similarity::ResultBlock::ResultBlock(int index, qint64 start):
   EAbstractAnalytic::Block(index),
   _start(start)
{
}






/*!
 * Append a pair to the result block's list of pairs.
 *
 * @param pair
 */
void Similarity::ResultBlock::append(const Pair& pair)
{
   _pairs.append(pair);
}






/*!
 * Write this block's data to the given data stream.
 *
 * @param stream
 */
void Similarity::ResultBlock::write(QDataStream& stream) const
{
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
