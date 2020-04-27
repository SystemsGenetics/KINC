#ifndef CLUSTER_FILTER_WORKBLOCK_H
#define CLUSTER_FILTER_WORKBLOCK_H
#include "corrpower.h"



/*!
 * This class implements the work block of the cluster_filter analytic.
 */
class CorrPowerFilter::WorkBlock : public EAbstractAnalyticBlock
{
    Q_OBJECT
public:
    /*!
     * Construct a new work block in an uninitialized null state.
     */
    explicit WorkBlock() = default;
    explicit WorkBlock(int index, qint64 start, qint64 size, Pairwise::Index startIndex);

    Pairwise::Index startIndex() const { return _startIndex; }
    qint64 size() const { return _size; }
    qint64 start() const {return _start;}

protected:
    virtual void write(QDataStream& stream) const override final;
    virtual void read(QDataStream& stream) override final;
private:
    /*!
     * The the number of the starting pair.
     */
    qint64 _start;
    /*!
     * The number of pairs to process.
     */
    qint64 _size;
    /*!
     * The starting index.
     */
    Pairwise::Index _startIndex;
};



#endif
