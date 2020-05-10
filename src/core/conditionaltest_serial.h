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
    int test(CCMatrix::Pair& ccmPair, qint32 clusterIndex, qint32& testIndex, qint32 featureIndex, qint32 labelIndex, QVector<QVector<double>>& pValues, QVector<QVector<double>>& r2);
    int prepAnxData(QString testLabel, int dataIndex, TestType testType);
    bool isEmpty(QVector<QVector<double>>& matrix);
    int clusterInfo(CCMatrix::Pair& ccmPair, int clusterIndex, QString label, TestType testType);

    // Binomial Tests
    double binomial();
    double testOne();
    double testTwo();

    // Hypergeometrix Test.
    double hypergeom(CCMatrix::Pair& ccmPair, int clusterIndex, QString test_label);

    // Regression Test
    void regression(QVector<QString> &amxInfo, CCMatrix::Pair& ccmPair, int clusterIndex, TestType testType, QVector<double>& results);
    double fTest(double chisq, gsl_matrix* X, gsl_vector* Y, gsl_matrix* cov, gsl_vector* C);

private:
    /*!
     * Pointer to the serials objects parent KNNAnalytic.
     */
    ConditionalTest* _base;
    /*!
     * Annotation matrix data for testing.
     */
    QVector<QString> _amxData;
    /*!
     * Category count.
     */
    qint32 _catCount {0};
    /*!
     * Category count in cluster.
     */
    qint32 _catInCluster {0};
    /*!
     * Size of the cluster.
     */
    qint32 _clusterSize {0};
};



#endif
