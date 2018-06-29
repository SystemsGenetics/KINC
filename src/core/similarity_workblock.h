#ifndef SIMILARITY_WORKBLOCK_H
#define SIMILARITY_WORKBLOCK_H
#include "similarity.h"



class Similarity::WorkBlock : public EAbstractAnalytic::Block
{
   Q_OBJECT
public:
   explicit WorkBlock() = default;
   explicit WorkBlock(int index, qint64 start, qint64 size);
   qint64 start() const { return _start; }
   qint64 size() const { return _size; }
protected:
   virtual void write(QDataStream& stream) const override final;
   virtual void read(QDataStream& stream) override final;
private:
   qint64 _start;
   qint64 _size;
};



#endif
