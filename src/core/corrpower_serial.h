#ifndef CLUSTER_FILTER_SERIAL_H
#define CLUSTER_FILTER_SERIAL_H
#include "corrpower.h"
#include "pairwise_clusteringmodel.h"
#include "pairwise_correlationmodel.h"



/*!
 * This class implements the serial working class of the similarity analytic.
 */
class CorrPowerFilter::Serial : public EAbstractAnalyticSerial
{
    Q_OBJECT
public:
    explicit Serial(CorrPowerFilter* parent);
    virtual std::unique_ptr<EAbstractAnalyticBlock> execute(const EAbstractAnalyticBlock* block) override final;
private:
     double pwr_r_test(double r, int n, double sig_level);
private:
    /*!
     * Pointer to the base analytic for this object.
     */
    CorrPowerFilter* _base;
};



#endif
