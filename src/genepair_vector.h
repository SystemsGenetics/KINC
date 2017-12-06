#ifndef GENEPAIR_VECTOR_H
#define GENEPAIR_VECTOR_H
#include <ace/core/AceCore.h>



namespace GenePair
{
   class Vector
   {
   public:
      Vector(qint64 geneX = 0, qint64 geneY = 0, qint64 cluster = 0);
      Vector(const Vector&) = default;
      Vector(Vector&&) = default;
      qint64 indent() const  { return (_geneX*(_geneX - 1)/2 + _geneY)*_maxClusterSize + _cluster; }
      qint64 geneX() const { return _geneX; }
      qint64 geneY() const { return _geneY; }
      qint64 cluster() const { return _cluster; }
      Vector& operator=(const Vector&) = default;
      Vector& operator=(Vector&&) = default;
      void operator++();
      bool operator==(const Vector& object)
         { return _geneX == object._geneX && _geneY == object._geneY
                  && _cluster == object._cluster; }
      bool operator!=(const Vector& object)
         { return _geneX != object._geneX || _geneY != object._geneY
                  || _cluster != object._cluster; }
      constexpr static qint64 _maxClusterSize {128};
   private:
      qint64 _geneX {0};
      qint64 _geneY {0};
      qint64 _cluster {0};
   };
}



#endif
