#ifndef GENEPAIR_PEARSON_H
#define GENEPAIR_PEARSON_H
#include "genepair_correlation.h"

namespace GenePair
{
   class Pearson : public Correlation
   {
   public:
      void initialize(ExpressionMatrix* /*input*/) {}
      QString getName() const { return "pearson"; }

   protected:
      float computeCluster(
         const QVector<Vector2>& data,
         const QVector<qint8>& labels,
         qint8 cluster,
         int minSamples
      );
   };
}

#endif
