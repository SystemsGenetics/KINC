#ifndef PAIRWISE_INDEX_H
#define PAIRWISE_INDEX_H
#include <ace/core/core.h>



namespace Pairwise
{
    /*!
     * This class implements the pairwise index, which provides a way to order
     * elements in a pairwise matrix and iterate through them. The pairwise index
     * uses row-major order and uses only the lower triangle of a matrix; that is,
     * it assumes that the row index is always greater than the column index.
     * Additionally, the pairwise index provides an "indent" value which can be
     * used to rank pairs that also have a cluster index; this value requires a
     * fixed upper bound on the number of clusters, which depends on the data
     * objects that use this class.
     */
    class Index
    {
    public:
        Index() = default;
        Index(qint32 x, qint32 y);
        Index(qint64 index);
        Index(const Index&) = default;
        Index(Index&&) = default;
        qint64 indent(qint8 cluster) const;
        qint32 getX() const { return _x; }
        qint32 getY() const { return _y; }
        Index& operator=(const Index&) = default;
        Index& operator=(Index&&) = default;
        void operator++();
        bool operator==(const Index& object) const
            { return _x == object._x && _y == object._y; }
        bool operator!=(const Index& object) const
            { return !(*this == object); }
        bool operator<(const Index& object) const
            { return _x < object._x || (_x == object._x && _y < object._y); }
        bool operator<=(const Index& object) const
            { return *this < object || *this == object; }
        bool operator>(const Index& object) const
            { return !(*this <= object); }
        bool operator>=(const Index& object) const
            { return !(*this < object); }
        /*!
         * The maximum number of clusters used to compute the indent value
         * of a pairwise index. Data objects which use the pairwise index should
         * never attempt to store more than this number of clusters in a single
         * pair.
         */
        constexpr static qint8 MAX_CLUSTER_SIZE {64};
    private:
        /*!
         * The row index.
         */
        qint32 _x {1};
        /*!
         * The column index.
         */
        qint32 _y {0};
    };
}



#endif
