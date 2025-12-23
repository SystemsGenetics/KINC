#ifndef IMPORTCSM_RESULTBLOCK_H
#define IMPORTCSM_RESULTBLOCK_H
#include "conditionaltest.h"



class ConditionalTest::ResultBlock : public EAbstractAnalyticBlock
{
    Q_OBJECT
public:
    explicit ResultBlock() = default;
    explicit ResultBlock(int index);
    const QVector<Pair>& pairs() const { return _pairs; }
    QVector<Pair>& pairs() { return _pairs; }
    void append(const Pair& pair);
protected:
    virtual void write(QDataStream& stream) const override final;
    virtual void read(QDataStream& stream) override final;
private:
    /*!
     * Row of dynamically populated cluster pvalues.
     */
    QVector<Pair> _pairs;
};



#endif
