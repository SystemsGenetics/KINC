#include "similarity_resultblock.h"






Similarity::ResultBlock::ResultBlock(int index, qint64 start):
   EAbstractAnalytic::Block(index),
   _start(start)
{
}






void Similarity::ResultBlock::append(const Pair& pair)
{
   _pairs.append(pair);
}






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
