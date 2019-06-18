#ifndef CLUSTER_FILTER_RESULTBLOCK_H
#define CLUSTER_FILTER_RESULTBLOCK_H
#include "corrpower.h"



/*!
 * This class implements the result block of the cluster_filter analytic.
 */
class CorrPowerFilter::ResultBlock : public EAbstractAnalyticBlock
{
   Q_OBJECT
public:
   /*!
    * Construct a new result block in an uninitialized null state.
    */
   explicit ResultBlock() = default;
   explicit ResultBlock(int index, qint64 start);
   template<class T> static QVector<T> makeVector(const T* data, int size);
   qint64 start() const { return _start; }
   const QVector<CPPair>& pairs() const { return _pairs; }
   QVector<CPPair>& pairs() { return _pairs; }
   void append(const CPPair& pair);
protected:
   virtual void write(QDataStream& stream) const override final;
   virtual void read(QDataStream& stream) override final;
private:
   /*!
    * The pairwise index of the first pair in the result block.
    */
   qint64 _start;
   /*!
    * The list of pairs that were processed.
    */
   QVector<CPPair> _pairs;
};



#endif
