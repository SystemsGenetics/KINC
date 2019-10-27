#ifndef CSM_PAIR_H
#define CSM_PAIR_H
#include <ace/core/core.h>
#include "conditionspecificclustersmatrix.h"
#include "pairwise_matrix_pair.h"



class CSMatrix::Pair : public Pairwise::Matrix::Pair
{
public:
    Pair(CSMatrix* matrix);
    Pair(const CSMatrix* matrix);
    Pair() = default;
    virtual ~Pair() = default;
public:
    //”cluster is our correlation stats”
    virtual void clearClusters() const;
    virtual void addCluster(int amount = 1) const;
    virtual int clusterSize() const ;
    virtual bool isEmpty() const;
    QString toString() const;
    const double& at(int cluster, int gene, QString type) const;
    double& at(int cluster, int gene, QString type);
    void addCluster(int amount, int size) const;
private:
    virtual void writeCluster(EDataStream& stream, int cluster);
    virtual void readCluster(const EDataStream& stream, int cluster) const;

    // Vectors for storing the p-values and R^2 values.
    // R^2 values are only for linear regression.
    mutable QVector<QVector<double>> _pValues;
    mutable QVector<QVector<double>> _r2;
    const CSMatrix* _cMatrix;
};



#endif
