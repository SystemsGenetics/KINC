#ifndef CORRELATIONMATRIX_PAIR_H
#define CORRELATIONMATRIX_PAIR_H
#include "correlationmatrix.h"
#include "pairwise_matrix_pair.h"



class CorrelationMatrix::Pair : public Pairwise::Matrix::Pair
{
public:
   Pair(CorrelationMatrix* matrix):
      Matrix::Pair(matrix),
      _cMatrix(matrix)
      {}
   Pair(const CorrelationMatrix* matrix):
      Matrix::Pair(matrix),
      _cMatrix(matrix)
      {}
   Pair() = default;
   virtual void clearClusters() const { _correlations.clear(); }
   virtual void addCluster(int amount = 1) const;
   virtual int clusterSize() const { return _correlations.size(); }
   virtual bool isEmpty() const { return _correlations.isEmpty(); }
   QString toString() const;
   const float& at(int cluster, int correlation) const
      { return _correlations.at(cluster).at(correlation); }
   float& at(int cluster, int correlation) { return _correlations[cluster][correlation]; }
private:
   virtual void writeCluster(EDataStream& stream, int cluster);
   virtual void readCluster(const EDataStream& stream, int cluster) const;
   mutable QVector<QVector<float>> _correlations;
   const CorrelationMatrix* _cMatrix;
};



#endif
