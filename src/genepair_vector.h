#ifndef GENEPAIR_VECTOR_H
#define GENEPAIR_VECTOR_H
#include <ace/core/AceCore.h>



namespace GenePair
{
   class Vector
   {
   public:
      Vector() = default;
      Vector(qint32 geneX, qint32 geneY);
      Vector(qint64 index);
      Vector(const Vector&) = default;
      Vector(Vector&&) = default;
      qint64 indent(qint8 cluster) const;
      qint32 geneX() const { return _geneX; }
      qint32 geneY() const { return _geneY; }
      Vector& operator=(const Vector&) = default;
      Vector& operator=(Vector&&) = default;
      void operator++();
      Vector operator++(int);
      bool operator==(const Vector& object)
         { return _geneX == object._geneX && _geneY == object._geneY; }
      bool operator!=(const Vector& object)
         { return !(*this == object); }
      bool operator<(const Vector& object)
         { return _geneX < object._geneX || (_geneX == object._geneX && _geneY < object._geneY); }
      bool operator<=(const Vector& object)
         { return *this < object || *this == object; }
      bool operator>(const Vector& object)
         { return !(*this <= object); }
      bool operator>=(const Vector& object)
         { return !(*this < object); }
      constexpr static qint8 MAX_CLUSTER_SIZE {64};
   private:
      qint32 _geneX {1};
      qint32 _geneY {0};
   };
}



#endif
