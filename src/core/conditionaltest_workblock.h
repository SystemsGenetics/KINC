#ifndef CSM_WORKBLOCK_H
#define CSM_WORKBLOCK_H
#include "conditionaltest.h"



class ConditionalTest::WorkBlock : public EAbstractAnalyticBlock
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
