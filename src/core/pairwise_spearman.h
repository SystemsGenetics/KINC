#ifndef PAIRWISE_SPEARMAN_H
#define PAIRWISE_SPEARMAN_H
#include "pairwise_correlationmodel.h"
#include "expressionmatrix.h"

namespace Pairwise
{
   /*!
    * This class implements the Spearman correlation model.
    */
   class Spearman : public CorrelationModel
   {
   public:
      static int nextPower2(int n);
      Spearman(ExpressionMatrix* emx);
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
      /*!
       * Workspace for the x data.
       */
      QVector<float> _x;
      /*!
       * Workspace for the y data.
       */
      QVector<float> _y;
      /*!
       * Workspace for the rank data.
       */
      QVector<int> _rank;
   };
}

#endif
