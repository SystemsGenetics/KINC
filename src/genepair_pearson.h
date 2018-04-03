#ifndef GENEPAIR_PEARSON_H
#define GENEPAIR_PEARSON_H
#include <ace/core/AceCore.h>

#include "genepair_correlation.h"

namespace GenePair
{
   class Pearson : public Correlation
   {
   public:
      double compute(
         ExpressionMatrix* input,
         Vector vector,
         const CCMatrix::Pair& pair, int cluster,
         int minSamples,
         int minExpression
      );
   
   private:
      QVector<double> _x;
      QVector<double> _y;
   };
}

#endif
