#ifndef GENEPAIR_CORRELATION_H
#define GENEPAIR_CORRELATION_H
#include <ace/core/AceCore.h>

#include "correlationmatrix.h"
#include "expressionmatrix.h"
#include "genepair_linalg.h"

namespace GenePair
{
   class Correlation
   {
   public:
      virtual void initialize(ExpressionMatrix* input) = 0;
      virtual QString getName() const = 0;

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
