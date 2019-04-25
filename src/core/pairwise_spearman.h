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
   public:
      Spearman(ExpressionMatrix* emx);
   protected:
      virtual float computeCluster(
         const float *x,
         const float *y,
         const QVector<qint8>& labels,
         qint8 cluster,
         int minSamples
      ) override final;
   private:
      template<typename T> void siftDown(QVector<float>& array, QVector<T>& extra, int start, int end);
      template<typename T> void heapSort(QVector<float>& array, QVector<T>& extra, int n);
   private:
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
