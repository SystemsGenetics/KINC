#ifndef CORRELATIONMATRIX_PAIR_H
#define CORRELATIONMATRIX_PAIR_H
#include "correlationmatrix.h"
#include "pairwise_matrix_pair.h"



/*!
 * This class implements the pairwise iterator for the correlation matrix data
 * object. This class extends the behavior of the base pairwise iterator to read
 * and write correlations.
 */
class CorrelationMatrix::Pair : public Pairwise::Matrix::Pair
{
public:
    Pair(CorrelationMatrix* matrix):
        Matrix::Pair(matrix),
        _cMatrix(matrix)
        {}
    Pair(const CorrelationMatrix* matrix):
        Matrix::Pair(matrix),
        _cMatrix(matrix)
        {}
    Pair() = default;
    virtual ~Pair() = default;
public:
    virtual void clearClusters() const { _correlations.clear(); }
    virtual void addCluster(int amount = 1) const;
    virtual int clusterSize() const { return _correlations.size(); }
    virtual bool isEmpty() const { return _correlations.isEmpty(); }
    QString toString() const;
    const float& at(int cluster) const { return _correlations.at(cluster); }
    float& at(int cluster) { return _correlations[cluster]; }
    const QVector<float>& correlations() { return _correlations; }
private:
    virtual void writeCluster(EDataStream& stream, int cluster);
    virtual void readCluster(const EDataStream& stream, int cluster) const;
    /*!
     * Array of correlations for the current pair.
     */
    mutable QVector<float> _correlations;
    /*!
     * Constant pointer to parent correlation matrix.
     */
    const CorrelationMatrix* _cMatrix;
};



#endif
