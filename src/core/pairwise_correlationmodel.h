#ifndef PAIRWISE_CORRELATIONMODEL_H
#define PAIRWISE_CORRELATIONMODEL_H
#include <ace/core/core.h>

#include "pairwise_index.h"

namespace Pairwise
{
   /*!
    * This class implements the abstract pairwise correlation model, which
    * takes a pairwise data array (with cluster labels) and computes a correlation
    * for each cluster in the data. The correlation metric must be implemented by
    * the inheriting class.
    */
   class CorrelationModel
   {
   public:
      ~CorrelationModel() = default;
   public:
      QVector<float> compute(
         const std::vector<float>& expressions,
         const Index& index,
         int K,
         const QVector<qint8>& labels,
         int minSamples
      );
   protected:
      virtual float computeCluster(
         const float *x,
         const float *y,
         const QVector<qint8>& labels,
         qint8 cluster,
         int minSamples
      ) = 0;
   };
}

#endif
