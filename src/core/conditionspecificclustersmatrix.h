#ifndef CSM_H
#define CSM_H
#include "pairwise_matrix.h"



class CSMatrix : public Pairwise::Matrix
{
    Q_OBJECT
public:
    class Pair;
public:
    virtual QAbstractTableModel* model() override final;
public:
    void initialize(const EMetaArray& geneNames, int maxClusterSize, int subheader);
    void initialize(const EMetaArray& features, const QVector<EMetaArray>& featureInfo, const QVector<EMetaArray>& data, int& numTests, QString fileName);

    int sampleSize() const;
    void setTestCount(qint32 newData);

    QString getTestName(int index) const;

    qint32 getTestCount();

private:
    class Model;
private:
    virtual void writeHeader() override final;
    virtual void readHeader() override final;
    qint32 _testcount {0};
    qint64 _sampleSize {0};
    constexpr static qint16 SUBHEADER_SIZE {12};
    Model* _model {nullptr};
};



#endif
