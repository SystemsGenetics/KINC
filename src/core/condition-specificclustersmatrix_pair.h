#ifndef CSCM_PAIR_H
#define CSCM_PAIR_H
#include <ace/core/core.h>
#include "condition-specificclustersmatrix.h"
#include "pairwise_matrix_pair.h"
//

class CSCM::Pair : public Pairwise::Matrix::Pair
{
public:
    Pair(CSCM* matrix);
    Pair(const CSCM* matrix);
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
    const CSCM* _cMatrix;
};

#endif // CSCM_PAIR_H
