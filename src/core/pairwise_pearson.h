#ifndef PAIRWISE_PEARSON_H
#define PAIRWISE_PEARSON_H
#include "pairwise_correlation.h"

namespace Pairwise
{
   /*!
    * This class implements the Pearson correlation model.
    */
   class Pearson : public Correlation
   {
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
