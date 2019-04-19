#ifndef SIMILARITY_WORKBLOCK_H
#define SIMILARITY_WORKBLOCK_H
#include "similarity.h"



/*!
 * This class implements the work block of the similarity analytic.
 */
class Similarity::WorkBlock : public EAbstractAnalyticBlock
{
   Q_OBJECT
public:
   /*!
    * Construct a new work block in an uninitialized null state.
    */
   explicit WorkBlock() = default;
   explicit WorkBlock(int index, qint64 start, qint64 size);
   qint64 start() const { return _start; }
   qint64 size() const { return _size; }
protected:
   virtual void write(QDataStream& stream) const override final;
   virtual void read(QDataStream& stream) override final;
private:
   /*!
    * The pairwise index of the first pair to process.
    */
   qint64 _start;
   /*!
    * The number of pairs to process.
    */
   qint64 _size;
};



#endif
