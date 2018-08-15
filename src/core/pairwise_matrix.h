#ifndef PAIRWISE_MATRIX_H
#define PAIRWISE_MATRIX_H
#include <ace/core/core.h>

#include "pairwise_index.h"



namespace Pairwise
{
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
      EMetadata geneNames() const;
   protected:
      virtual void writeHeader() = 0;
      virtual void readHeader() = 0;
      void initialize(const EMetadata& geneNames, int maxClusterSize, int dataSize, int offset);
   private:
      void write(Index index, qint8 cluster);
      Index getPair(qint64 index, qint8* cluster) const;
      qint64 findPair(qint64 indent, qint64 first, qint64 last) const;
      void seekPair(qint64 index) const;
      constexpr static int _headerSize {30};
      constexpr static int _itemHeaderSize {9};
      qint32 _geneSize {0};
      qint32 _maxClusterSize {0};
      qint32 _dataSize {0};
      qint64 _pairSize {0};
      qint64 _clusterSize {0};
      qint16 _offset {0};
      qint64 _lastWrite {-2};
   };
}



#endif
