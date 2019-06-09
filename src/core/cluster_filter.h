#ifndef CLUSTER_FILTER_H
#define CLUSTER_FILTER_H

#include <ace/core/core.h>

#include "ccmatrix_pair.h"
#include "ccmatrix.h"
#include "correlationmatrix_pair.h"
#include "correlationmatrix.h"
#include "expressionmatrix.h"

class ClusterFilter : public EAbstractAnalytic
{
    Q_OBJECT
 public:
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
 private:
    /**
     * Workspace variables to write to the output file
     */
    QTextStream _stream;
    CCMatrix::Pair _ccmPair;
    CorrelationMatrix::Pair _cmxPair;
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
     * Pointer to the output text file.
     */
    QFile* _output {nullptr};
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
};

#endif
