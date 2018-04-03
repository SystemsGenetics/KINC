#ifndef GENEPAIR_CORRELATION_H
#define GENEPAIR_CORRELATION_H
#include <ace/core/AceCore.h>

#include "ccmatrix.h"
#include "expressionmatrix.h"
#include "genepair_vector.h"

namespace GenePair
{
   class Correlation
   {
   public:
      virtual double compute(
         ExpressionMatrix* input,
         Vector vector,
         const CCMatrix::Pair& pair, int cluster,
         int minSamples,
         int minExpression
      ) = 0;

   protected:
      void fetchData(
         ExpressionMatrix* input,
         Vector vector,
         const CCMatrix::Pair& pair, int cluster,
         int minExpression,
         QVector<double>& x,
         QVector<double>& y
      );
   };
}

#endif
