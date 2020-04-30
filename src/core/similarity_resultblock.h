#ifndef SIMILARITY_RESULTBLOCK_H
#define SIMILARITY_RESULTBLOCK_H
#include "similarity.h"



/*!
 * This class implements the result block of the similarity analytic.
 */
class Similarity::ResultBlock : public EAbstractAnalyticBlock
{
    Q_OBJECT
public:
    /*!
     * Construct a new result block in an uninitialized null state.
     */
    explicit ResultBlock() = default;
    explicit ResultBlock(int index);
    template<class T> static QVector<T> makeVector(const T* data, int size);
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



/*!
 * Create a vector from the given pointer and size. The contents of the
 * pointer are copied into the vector.
 *
 * @param data
 * @param size
 */
template<class T>
QVector<T> Similarity::ResultBlock::makeVector(const T* data, int size)
{
    QVector<T> v(size);

    memcpy(v.data(), data, size * sizeof(T));
    return v;
}



#endif
