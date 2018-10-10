#ifndef PAIRWISE_MATRIX_PAIR_H
#define PAIRWISE_MATRIX_PAIR_H
#include "pairwise_matrix.h"



namespace Pairwise
{
   /*!
    * This class implements the pairwise iterator for the pairwise matrix
    * data object. The pairwise iterator can read from or write to any pair in
    * the pairwise matrix, or it can simply iterate through each pair. The
    * iterator stores only one pair in memory at a time.
    */
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
      void write(const Index& index);
      void read(const Index& index) const;
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
      /*!
       * Pointer to the parent pairwise matrix.
       */
      Matrix* _matrix {nullptr};
      /*!
       * Constant pointer to the parent pairwise matrix.
       */
      const Matrix* _cMatrix;
      /*!
       * The iterator's current position in the pairwise matrix.
       */
      mutable qint64 _rawIndex {0};
      /*!
       * Pairwise index corresponding to the iterator's position.
       */
      mutable Index _index;
   };
}



#endif
