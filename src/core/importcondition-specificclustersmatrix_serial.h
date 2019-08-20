#ifndef IMPORTCSCM_SERIAL_H
#define IMPORTCSCM_SERIAL_H
#include <ace/core/core.h>
#include "importcondition-specificclustersmatrix.h"

class importCSCM::Serial : public EAbstractAnalyticSerial
{
   Q_OBJECT
public:
   explicit Serial(importCSCM* parent);
   virtual std::unique_ptr<EAbstractAnalyticBlock> execute(const EAbstractAnalyticBlock* block) override final;

    //helper functions
    int test(CorrelationMatrix::Pair cmxPair, CCMatrix::Pair& ccmPair, qint32 clusterIndex, qint32& testIndex, qint32 featureIndex, qint32 labelIndex, QVector<QVector<double>>& pValues);
    int prepAnxData(QString testLabel, int dataIndex);
    bool isEmpty(QVector<QVector<double>>& matrix);
    int clusterInfo(CCMatrix::Pair& ccmPair, int clusterIndex, QString label);

    //Binomial Tests
    double binomial(double alpha);
    double testOne();
    double testTwo();

    //Regresion Test
    double regresion(QVector<QString>& anxInfo, CCMatrix::Pair& ccmPair, int clusterIndex);

private:
   /*!
   *  Pointer to the serials objects parent KNNAnalytic.
   */
   importCSCM* _base;
   /*!
   *  Annotation matrix data for testing.
   */
   QVector<QString> _anxData;
   /*!
   *  Catagory count.
   */
   qint32 _catCount {0};
   /*!
   *  Catagory count in cluster.
   */
   qint32 _catInCount {0};
   /*!
   *  Catagory count in cluster.
   */
   qint32 _clusterInMask {0};
};

#endif // IMPORTCSCM_SERIAL_H
