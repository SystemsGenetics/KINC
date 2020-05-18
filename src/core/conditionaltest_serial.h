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
    void hypergeom(QVector<QString> amx_column, CCMatrix::Pair& ccmPair, int clusterIndex,
                   int featureIndex, int labelIndex,  double& results);
    void regression(QVector<QString> amx_column, CCMatrix::Pair& ccmPair, int clusterIndex,
                    int featureIndex, QVector<double>& results);

private:
    /*!
     * Pointer to the base analytic for this object.
     */
    ConditionalTest* _base;
};



#endif
