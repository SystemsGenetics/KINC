#ifndef PAIRWISE_MATRIX_H
#define PAIRWISE_MATRIX_H
#include <ace/core/core.h>

#include "pairwise_index.h"



namespace Pairwise
{
   /*!
    * This class implements the abstract pairwise matrix data object, which can
    * be extended to represent any pairwise matrix. Both the rows and columns
    * correspond to genes, and each element (i, j) in the matrix contains
    * pairwise data for genes i and j. This pairwise data can have multiple clusters,
    * and the structure of a "pair-cluster" is defined by the inheriting class.
    * This class stores matrix data as an ordered list of indexed pairs; therefore,
    * pairwise data must be written in order and it should be sparse for the
    * storage format to be efficient.
    */
   class Matrix : public EAbstractData
   {
   public:
      class Pair;
   public:
      virtual qint64 dataEnd() const override final;
      virtual void readData() override final;
      virtual void writeNewData() override final;
      virtual void finish() override final;
   public:
      int geneSize() const { return _geneSize; }
      int maxClusterSize() const { return _maxClusterSize; }
      qint64 size() const { return _pairSize; }
      EMetaArray geneNames() const;
   protected:
      virtual void writeHeader() = 0;
      virtual void readHeader() = 0;
      void initialize(const EMetaArray& geneNames, int maxClusterSize, int dataSize, int offset);
   private:
      void write(const Index& index, qint8 cluster);
      Index getPair(qint64 index, qint8* cluster) const;
      qint64 findPair(qint64 indent, qint64 first, qint64 last) const;
      void seekPair(qint64 index) const;
      /*!
       * The size (in bytes) of the header at the beginning of the file. The header
       * consists of the gene size, max cluster size, pairwise data size, total
       * number of pairs, total number of clusters, and sub-header offset.
       */
      constexpr static int _headerSize {30};
      /*!
       * The size (in bytes) of the pairwise header. The item header size consists
       * of the row and column index of the pair.
       */
      constexpr static int _itemHeaderSize {9};
      /*!
       * The number of genes in the pairwise matrix.
       */
      qint32 _geneSize {0};
      /*!
       * The maximum number of clusters allowed for each pair in the matrix.
       */
      qint32 _maxClusterSize {0};
      /*!
       * The size (in bytes) of a pairwise data element.
       */
      qint32 _dataSize {0};
      /*!
       * The total number of pairs in the matrix.
       */
      qint64 _pairSize {0};
      /*!
       * The total number of clusters (across all pairs) in the matrix.
       */
      qint64 _clusterSize {0};
      /*!
       * The size (in bytes) of the sub-header, which occurs after the header
       * and can be used by an inheriting class.
       */
      qint16 _subHeaderSize {0};
      /*!
       * The index of the last pair that was written to the matrix.
       */
      qint64 _lastWrite {-2};
   };
}



#endif
