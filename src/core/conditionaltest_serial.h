#ifndef ConditionalTest_SERIAL_H
#define ConditionalTest_SERIAL_H
#include <ace/core/core.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_randist.h>
#include "conditionaltest.h"



class ConditionalTest::Serial : public EAbstractAnalyticSerial
{
    Q_OBJECT
public:
    explicit Serial(ConditionalTest* parent);
    virtual std::unique_ptr<EAbstractAnalyticBlock> execute(const EAbstractAnalyticBlock* block) override final;

    // helper functions
    bool isEmpty(QVector<QVector<double>>& matrix);

    // Statistical Tests
    void hypergeom(
        const QVector<QString>& amx_column,
        const CCMatrix::Pair& ccmPair,
        int clusterIndex,
        int featureIndex,
        int labelIndex,
        double& results);
    
    void regression(
        const QVector<QString>& amx_column,
        const CCMatrix::Pair& ccmPair,
        int clusterIndex,
        int featureIndex,
        QVector<double>& results);

private:
    /*!
     * Pointer to the base analytic for this object.
     */
    ConditionalTest* _base;

    // Performs power analysis for multiple-linear regression.
    double pwr_f2_test(int u, int v, double f2, double sig_level);
};



#endif
