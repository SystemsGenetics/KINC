#ifndef CLUSTER_FILTER_SERIAL_H
#define CLUSTER_FILTER_SERIAL_H
#include "cluster_filter.h"
#include "pairwise_clusteringmodel.h"
#include "pairwise_correlationmodel.h"



/*!
 * This class implements the serial working class of the similarity analytic.
 */
class ClusterFilter::Serial : public EAbstractAnalyticSerial
{
   Q_OBJECT
public:
   explicit Serial(ClusterFilter* parent);
   virtual std::unique_ptr<EAbstractAnalyticBlock> execute(const EAbstractAnalyticBlock* block) override final;
private:

private:
   /*!
    * Pointer to the base analytic for this object.
    */
   ClusterFilter* _base;

};


#endif
