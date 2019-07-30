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
      void siftDown(QVector<float>& array, QVector<float>& extra, int start, int end);
      void heapSort(QVector<float>& array, QVector<float>& extra, int n);
      void computeRank(QVector<float>& array, int n);
   private:
      /*!
       * Workspace for the x rank data.
       */
      QVector<float> _x_rank;
      /*!
       * Workspace for the y rank data.
       */
      QVector<float> _y_rank;
   };
}

#endif
