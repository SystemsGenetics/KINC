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
    const double& at(int cluster, int gene) const;
    double& at(int cluster, int gene);
    void addCluster(int amount, int size) const;
private:
    virtual void writeCluster(EDataStream& stream, int cluster);
    virtual void readCluster(const EDataStream& stream, int cluster) const;

    mutable QVector<QVector<double>> _pValues;
    const CSMatrix* _cMatrix;
};



#endif
