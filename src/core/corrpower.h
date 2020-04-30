#ifndef CLUSTER_FILTER_H
#define CLUSTER_FILTER_H

#include <ace/core/core.h>

#include "ccmatrix_pair.h"
#include "ccmatrix.h"
#include "correlationmatrix_pair.h"
#include "correlationmatrix.h"
#include "expressionmatrix.h"
#include "pairwise_clusteringmodel.h"



/*!
 * This class implements the cluster filter analytic. This analytic takes a
 * correlation matrix and cluster matrix and filters out clusters that do not
 * pass specific tests. This analytic can use MPI.
 */
class CorrPowerFilter : public EAbstractAnalytic
{
    Q_OBJECT
public:
    /*!
     * Defines the pair structure used to send results in result blocks.
     */
    struct Pair
    {
        /*!
         * The pairwise index of the pair.
         */
        Pairwise::Index index;
        /*!
         * The cluster labels for a pair.
         */
        QVector<qint8> labels;
        /*!
         * The correlation for each cluster in a pair.
         */
        QVector<float> correlations;
        /*!
         * The cluster indexes to keep.
         */
        QVector<int> keep;
    };
    class Input;
    class WorkBlock;
    class ResultBlock;
    class Serial;
public:
    virtual int size() const override final;
    virtual std::unique_ptr<EAbstractAnalyticBlock> makeWork(int index) const override final;
    virtual std::unique_ptr<EAbstractAnalyticBlock> makeWork() const override final;
    virtual std::unique_ptr<EAbstractAnalyticBlock> makeResult() const override final;
    virtual void process(const EAbstractAnalyticBlock* result) override final;
    virtual EAbstractAnalyticInput* makeInput() override final;
    virtual EAbstractAnalyticSerial* makeSerial() override final;
    virtual void initialize() override final;
    virtual void initializeOutputs() override final;
private:
    /*!
     * Pointer to the input cluster matrix.
     */
    CCMatrix* _ccm {nullptr};
    /*!
     * Pointer to the input correlation matrix.
     */
    CorrelationMatrix* _cmx {nullptr};
    /*!
     * Pointer to the output cluster matrix.
     */
    CCMatrix* _ccmOut {nullptr};
    /*!
     * Pointer to the output correlation matrix.
     */
    CorrelationMatrix* _cmxOut {nullptr};
    /*!
     * The significance level (i.e. Type I error rate, alpha) for the power test.
     */
    double _powerThresholdAlpha {0.001};
    /*!
     * The power value (i.e. 1 minus Type II error rate, 1 minus beta) for the power test.
     */
    double _powerThresholdPower {0.8};
    /*!
     * The number of pairs to process in each work block.
     */
    int _workBlockSize {0};
};



#endif
