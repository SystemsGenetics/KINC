#ifndef GENEPAIR_BASE_H
#define GENEPAIR_BASE_H
#include <ace/core/AceCore.h>

#include "genepair_vector.h"



namespace GenePair
{
   class Base : public EAbstractData
   {
   public:
      class Pair
      {
      public:
         Pair(Base* matrix);
         Pair(const Base* matrix);
         virtual void clear() = 0;
         virtual void addCluster() = 0;
         void write(Vector index);
         void read(Vector index) const;
         void readFirst() const;
         void readNext() const;
         bool hasNext() const;
         int clusterSize() const;
         bool isEmpty() const;
         const Vector& index() const;
      protected:
         virtual void writeCluster(EDataStream& stream, int cluster) = 0;
         virtual void readCluster(EDataStream& stream, int cluster) = 0;
      private:
         Base* _matrix {nullptr};
         const Base* _cMatrix;
         Vector _index;
      };
      virtual void readData() override final;
      virtual quint64 getDataEnd() const override final
         { return _headerSize + _offset + _pairSize*(_dataSize + _itemHeaderSize); }
      virtual void newData() override final;
      virtual void prepare(bool) override final {}
      virtual void finish() override final { Base::newData(); }
      qint64 size() const { return _pairSize; }
   protected:
      virtual void writeHeader() = 0;
      virtual void readHeader() = 0;
      void initialize(int geneSize, int dataSize, int offset);
   private:
      void write(Vector index);
      Vector getPair(qint64 index) const;
      qint64 findPair(qint64 indent, qint64 first, qint64 last) const;
      void seekPair(qint64 index) const;
      constexpr static int _headerSize {24};
      constexpr static int _itemHeaderSize {9};
      qint32 _geneSize {0};
      qint32 _dataSize {0};
      qint64 _pairSize {0};
      qint64 _rawPairSize {0};
      qint16 _offset {0};
      qint64 _lastWrite;
   };
}



#endif
