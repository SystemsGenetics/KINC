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
 * This class implements the cluster filter analytic. This analytic takes an
 * expression matrix, correlation matrix and cluster composition matrix allows
 * the user to filter out clusters that do not pass specific tests.
 *
 * This analytic can use MPI.
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
        * The number of clusters in a pair.
        */
       qint8 K;
       /*!
        * The cluster labels for a pair.
        */
       QVector<qint8> labels;
       /*!
        * The correlation for each cluster in a pair.
        */
       QVector<float> correlations;
       /*!
        * The x/y coordinates in the CCM/CMX matricies that this pair belongs to.
        */
       qint32 x_index;
       qint32 y_index;

    };
    class Input;
    class WorkBlock;
    class ResultBlock;
    class Serial;
    virtual int size() const override final;
    virtual std::unique_ptr<EAbstractAnalyticBlock> makeWork(int index) const override final;
    virtual std::unique_ptr<EAbstractAnalyticBlock> makeWork() const override final;
    virtual std::unique_ptr<EAbstractAnalyticBlock> makeResult() const override final;
    virtual void process(const EAbstractAnalyticBlock* result) override final;
    virtual EAbstractAnalyticInput* makeInput() override final;
    virtual EAbstractAnalyticSerial* makeSerial() override final;
    virtual void initialize() override final;
    virtual void initializeOutputs() override final;
    static qint64 totalPairs(const ExpressionMatrix* emx);
 private:
    /*!
     * Pointer to the input expression matrix.
     */
    ExpressionMatrix* _emx {nullptr};
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
     * Whether to remove outliers before clustering.
     */
    bool _doCorrelationPowerFilter {false};
    /*!
     * The minimum (absolute) correlation threshold to save a correlation.
     */
    double _powerThresholdAlpha {0.001};
    /*!
     * The minimum (absolute) correlation threshold to save a correlation.
     */
    double _powerThresholdPower {0.8};
    /*!
     * The number of pairs to process in each work block.
     */
    int _workBlockSize {0};
};

#endif
