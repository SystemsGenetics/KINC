#ifndef PAIRWISE_CORRELATION_H
#define PAIRWISE_CORRELATION_H
#include <ace/core/core.h>

#include "correlationmatrix.h"
#include "expressionmatrix.h"
#include "pairwise_linalg.h"

namespace Pairwise
{
   class Correlation
   {
   public:
      QVector<float> compute(
         const QVector<Vector2>& data,
         int K,
         const QVector<qint8>& labels,
         int minSamples
      );

   protected:
      virtual float computeCluster(
         const QVector<Vector2>& data,
         const QVector<qint8>& labels,
         qint8 cluster,
         int minSamples
      ) = 0;
   };
}

#endif
