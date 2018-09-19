#ifndef PAIRWISE_SPEARMAN_H
#define PAIRWISE_SPEARMAN_H
#include "pairwise_correlation.h"

namespace Pairwise
{
   class Spearman : public Correlation
   {
   public:
      static int nextPower2(int n);

      Spearman(ExpressionMatrix* input);

   protected:
      float computeCluster(
         const QVector<Vector2>& data,
         const QVector<qint8>& labels,
         qint8 cluster,
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
