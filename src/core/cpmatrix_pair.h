#ifndef CPMATRIX_PAIR_H
#define CPMATRIX_PAIR_H
#include "cpmatrix.h"
#include "pairwise_matrix_pair.h"



/*!
 * This class implements the pairwise iterator for the correlation matrix data
 * object. This class extends the behavior of the base pairwise iterator to read
 * and write components.
 */
class CPMatrix::Pair : public Pairwise::Matrix::Pair
{
public:
   Pair(CPMatrix* matrix):
      Matrix::Pair(matrix),
      _cMatrix(matrix)
      {}
   Pair(const CPMatrix* matrix):
      Matrix::Pair(matrix),
      _cMatrix(matrix)
      {}
   Pair() = default;
   virtual ~Pair() = default;
public:
   virtual void clearClusters() const { _components.clear(); }
   virtual void addCluster(int amount = 1) const;
   virtual int clusterSize() const { return _components.size(); }
   virtual bool isEmpty() const { return _components.isEmpty(); }
   QString toString() const;
   const Component& at(int cluster) const { return _components.at(cluster); }
   Component& at(int cluster) { return _components[cluster]; }
   const QVector<Component>& components() { return _components; }
private:
   virtual void writeCluster(EDataStream& stream, int cluster);
   virtual void readCluster(const EDataStream& stream, int cluster) const;
   /*!
    * Array of components for the current pair.
    */
   mutable QVector<Component> _components;
   /*!
    * Constant pointer to parent correlation matrix.
    */
   const CPMatrix* _cMatrix;
};



#endif
