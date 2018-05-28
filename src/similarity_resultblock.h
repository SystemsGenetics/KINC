#ifndef SIMILARITY_RESULTBLOCK_H
#define SIMILARITY_RESULTBLOCK_H
#include "similarity.h"



class Similarity::ResultBlock : public EAbstractAnalytic::Block
{
   Q_OBJECT
public:
   explicit ResultBlock() = default;
   explicit ResultBlock(int index, qint64 start);
   qint64 start() const { return _start; }
   const QVector<Pair>& pairs() const { return _pairs; }
   QVector<Pair>& pairs() { return _pairs; }
   void append(const Pair& pair);
protected:
   virtual void write(QDataStream& stream) const override final;
   virtual void read(QDataStream& stream) override final;
private:
   qint64 _start;
   QVector<Pair> _pairs;
};



#endif
