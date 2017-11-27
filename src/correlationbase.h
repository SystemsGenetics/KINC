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
      bool operator<(const Iterator& object) const { return indent() < object.indent(); }
      bool operator>(const Iterator& object) const { return indent() > object.indent(); }
      bool operator<=(const Iterator& object) const { return indent() <= object.indent(); }
      bool operator>=(const Iterator& object) const { return indent() >= object.indent(); }
   private:
      int _x {1};
      int _y {0};
   };
protected:
   void initialize(int geneSize, int dataSize, int offset);
   void write(Iterator correlation);
   void read();
   bool seek(Iterator correlation) const;
private:
   int _geneSize {0};
   int _dataSize {0};
   int _offset {0};
   Iterator _lastWrite;
};



#endif
