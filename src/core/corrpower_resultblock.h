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
    explicit ResultBlock(int index);
    const QVector<Pair>& pairs() const { return _pairs; }
    QVector<Pair>& pairs() { return _pairs; }
    void append(const Pair& pair);
protected:
    virtual void write(QDataStream& stream) const override final;
    virtual void read(QDataStream& stream) override final;
private:
    /*!
     * The list of pairs that were processed.
     */
    QVector<Pair> _pairs;
};



#endif
