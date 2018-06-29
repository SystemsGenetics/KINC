#ifndef PAIRWISE_INDEX_H
#define PAIRWISE_INDEX_H
#include <ace/core/core.h>



namespace Pairwise
{
   class Index
   {
   public:
      Index() = default;
      Index(qint32 x, qint32 y);
      Index(qint64 index);
      Index(const Index&) = default;
      Index(Index&&) = default;
      qint64 indent(qint8 cluster) const;
      qint32 getX() const { return _x; }
      qint32 getY() const { return _y; }
      Index& operator=(const Index&) = default;
      Index& operator=(Index&&) = default;
      void operator++();
      Index operator++(int);
      bool operator==(const Index& object) const
         { return _x == object._x && _y == object._y; }
      bool operator!=(const Index& object)
         { return !(*this == object); }
      bool operator<(const Index& object)
         { return _x < object._x || (_x == object._x && _y < object._y); }
      bool operator<=(const Index& object)
         { return *this < object || *this == object; }
      bool operator>(const Index& object)
         { return !(*this <= object); }
      bool operator>=(const Index& object)
         { return !(*this < object); }
      constexpr static qint8 MAX_CLUSTER_SIZE {64};
   private:
      qint32 _x {1};
      qint32 _y {0};
   };
}



#endif
