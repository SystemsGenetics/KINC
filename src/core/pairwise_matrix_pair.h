#ifndef PAIRWISE_MATRIX_PAIR_H
#define PAIRWISE_MATRIX_PAIR_H
#include "pairwise_matrix.h"



namespace Pairwise
{
   class Matrix::Pair
   {
   public:
      Pair(Matrix* matrix):
         _matrix(matrix),
         _cMatrix(matrix),
         _index({matrix->_geneSize,0})
         {}
      Pair(const Matrix* matrix):
         _cMatrix(matrix),
         _index({matrix->_geneSize,0})
         {}
      Pair() = default;
      Pair(const Pair&) = default;
      Pair(Pair&&) = default;
      virtual void clearClusters() const = 0;
      virtual void addCluster(int amount = 1) const = 0;
      virtual int clusterSize() const = 0;
      virtual bool isEmpty() const = 0;
      void write(Index index);
      void read(Index index) const;
      void reset() const { _rawIndex = 0; };
      void readNext() const;
      bool hasNext() const { return _rawIndex != _cMatrix->_clusterSize; }
      const Index& index() const { return _index; }
      Pair& operator=(const Pair&) = default;
      Pair& operator=(Pair&&) = default;
   protected:
      virtual void writeCluster(EDataStream& stream, int cluster) = 0;
      virtual void readCluster(const EDataStream& stream, int cluster) const = 0;
   private:
      Matrix* _matrix {nullptr};
      const Matrix* _cMatrix;
      mutable qint64 _rawIndex {0};
      mutable Index _index;
   };
}



#endif
