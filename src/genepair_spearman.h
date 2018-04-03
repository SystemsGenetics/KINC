#ifndef GENEPAIR_SPEARMAN_H
#define GENEPAIR_SPEARMAN_H
#include <ace/core/AceCore.h>

#include "genepair_correlation.h"

namespace GenePair
{
   class Spearman : public Correlation
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
      QVector<double> _work;
   };
}

#endif
