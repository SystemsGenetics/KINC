#ifndef GENEPAIR_SPEARMAN_H
#define GENEPAIR_SPEARMAN_H
#include "genepair_correlation.h"

namespace GenePair
{
   class Spearman : public Correlation
   {
   public:
      void initialize(ExpressionMatrix* input, CorrelationMatrix* output);
      float compute(
         const QVector<Vector2>& data,
         const QVector<qint8>& labels, qint8 cluster,
         int minSamples
      );

   private:
      void bitonicSort(int size, QVector<float>& sortList, QVector<float>& extraList);
      void bitonicSort(int size, QVector<float>& sortList, QVector<int>& extraList);

      QVector<float> _x;
      QVector<float> _y;
      QVector<float> _rank;
   };
}

#endif
