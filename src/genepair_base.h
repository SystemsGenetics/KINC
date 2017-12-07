#ifndef GENEPAIR_BASE_H
#define GENEPAIR_BASE_H
#include <ace/core/AceCore.h>

#include "genepair_vector.h"



namespace GenePair
{
   class Base : public EAbstractData
   {
   public:
      virtual void readData() override;
      virtual quint64 getDataEnd() const override
         { return _headerSize + _offset + _pairSize*(_dataSize + _itemHeaderSize); }
      virtual void newData() override;
      virtual void prepare(bool preAllocate) override;
      virtual void finish() override { Base::newData(); }
   protected:
      void initialize(int geneSize, int dataSize, int offset);
      void write(Vector index);
      qint64 pairSize() const { return _pairSize; }
      Vector getPair(qint64 index) const;
      qint64 findPair(Vector index) const { return findPair(index.indent(),0,_pairSize - 1); }
   private:
      qint64 findPair(qint64 indent, qint64 first, qint64 last) const;
      void seekPair(qint64 index) const;
      constexpr static int _headerSize {18};
      constexpr static int _itemHeaderSize {16};
      qint32 _geneSize {0};
      qint32 _dataSize {0};
      qint64 _pairSize {0};
      qint16 _offset {0};
      Vector _lastWrite;
   };
}



#endif
