#ifndef CSM_WORKBLOCK_H
#define CSM_WORKBLOCK_H
#include <ace/core/core.h>
#include "conditionaltest.h"



class ConditionalTest::WorkBlock : public EAbstractAnalyticBlock
{
    Q_OBJECT
public:
    /*!
     * Creates an uninitialized work block
     */
    explicit WorkBlock() = default;
    explicit WorkBlock(int index, Pairwise::Index start, qint64 startpair, qint64 size);

    Pairwise::Index start() const { return _start; }
    qint64 size() const { return _size; }
    qint64 startpair() const {return _startpair;}
protected:
    virtual void write(QDataStream& stream) const override final;
    virtual void read(QDataStream& stream) override final;
private:
    /*!
     * The pairwise index of the first pair to process.
     */
    Pairwise::Index _start;
    /*!
     * The the number of the starting pair.
     */
    qint64 _startpair;
    /*!
     * The number of pairs to process.
     */
    qint64 _size;
};



#endif
