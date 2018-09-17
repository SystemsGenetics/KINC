#ifndef CCMATRIX_PAIR_H
#define CCMATRIX_PAIR_H
#include "ccmatrix.h"
#include "pairwise_matrix_pair.h"



/*!
 * This class implements the pairwise iterator for the cluster matrix data
 * object. This class extends the behavior of the base pairwise iterator to read
 * and write sample masks.
 */
class CCMatrix::Pair : public Pairwise::Matrix::Pair
{
public:
   Pair(CCMatrix* matrix):
      Matrix::Pair(matrix),
      _cMatrix(matrix)
      {}
   Pair(const CCMatrix* matrix):
      Matrix::Pair(matrix),
      _cMatrix(matrix)
      {}
   Pair() = default;
   virtual void clearClusters() const { _sampleMasks.clear(); }
   virtual void addCluster(int amount = 1) const;
   virtual int clusterSize() const { return _sampleMasks.size(); }
   virtual bool isEmpty() const { return _sampleMasks.isEmpty(); }
   QString toString() const;
   const qint8& at(int cluster, int sample) const { return _sampleMasks.at(cluster).at(sample); }
   qint8& at(int cluster, int sample) { return _sampleMasks[cluster][sample]; }
private:
   virtual void writeCluster(EDataStream& stream, int cluster);
   virtual void readCluster(const EDataStream& stream, int cluster) const;
   /*!
    * Array of sample masks for the current pair.
    */
   mutable QVector<QVector<qint8>> _sampleMasks;
   /*!
    * Constant pointer to parent cluster matrix.
    */
   const CCMatrix* _cMatrix;
};



#endif
