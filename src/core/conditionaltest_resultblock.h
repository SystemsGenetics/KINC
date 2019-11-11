#ifndef IMPORTCSM_RESULTBLOCK_H
#define IMPORTCSM_RESULTBLOCK_H
#include <ace/core/core.h>
#include "conditionaltest.h"



class ConditionalTest::ResultBlock : public EAbstractAnalyticBlock
{
   Q_OBJECT
public:
    explicit ResultBlock() = default;
    explicit ResultBlock(int index);
    explicit ResultBlock(int index, int numTests, qint64 start);

    qint64 start() const { return _start; }
    const QVector<CSMPair>& pairs() const { return _pairs; }
    QVector<CSMPair>& pairs() { return _pairs; }
    void append(const CSMPair& pair);

protected:
    virtual void write(QDataStream& stream) const override final;
    virtual void read(QDataStream& stream) override final;
private:
    /*!
    *  Number of each test to be conducted on each cluster.
    */
    int _numTests{0};
    /*!
     * The pairwise index of the first pair in the result block.
     */
    qint64 _start{0};
    /*!
    *  Row of dynamically populated cluster pvalues.
    */
    QVector<CSMPair> _pairs;
};



#endif
