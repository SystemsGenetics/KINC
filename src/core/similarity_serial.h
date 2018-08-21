#ifndef SIMILARITY_SERIAL_H
#define SIMILARITY_SERIAL_H
#include "similarity.h"



class Similarity::Serial : public EAbstractAnalytic::Serial
{
   Q_OBJECT
public:
   explicit Serial(Similarity* parent);
   virtual std::unique_ptr<EAbstractAnalytic::Block> execute(const EAbstractAnalytic::Block* block) override final;
private:
   int fetchPair(Pairwise::Index index, QVector<Pairwise::Vector2>& data, QVector<qint8>& labels);
   void markOutliers(const QVector<Pairwise::Vector2>& data, int N, int j, QVector<qint8>& labels, qint8 cluster, qint8 marker);

   Similarity* _base;
};



#endif
