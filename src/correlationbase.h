#ifndef CORRELATIONBASE_H
#define CORRELATIONBASE_H
#include <ace/core/AceCore.h>



class CorrelationBase : public EAbstractData
{
public:
   class Iterator
   {
   public:
      Iterator() = default;
      Iterator(int x, int y);
      int x() const { return _x; }
      int y() const { return _y; }
      qint64 indent() const { return (qint64)_x*((qint64)_x - 1)/2 + (qint64)_y; }
      void operator--();
      void operator++();
      bool operator!=(const Iterator& object) const { return _x != object._x || _y != object._y; }
      bool operator==(const Iterator& object) const { return _x == object._x && _y == object._y; }
   private:
      int _x {1};
      int _y {0};
   };
   void finish() override final;
protected:
   void initialize(int geneSize, int dataSize, int offset);
   void write(Iterator index);
   void read();
   bool findCorrelation(Iterator index) const
      { return findCorrelation(index.indent(),0,_correlationSize - 1); }
private:
   bool findCorrelation(qint64 indent, int first, int last) const;
   void seekCorrelation(int index) const;
   constexpr static int _headerSize {18};
   qint32 _geneSize {0};
   qint32 _dataSize {0};
   qint64 _correlationSize {0};
   qint16 _offset {0};
   Iterator _lastWrite;
};



#endif
