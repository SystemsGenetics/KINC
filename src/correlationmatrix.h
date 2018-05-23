#ifndef CORRELATIONMATRIX_H
#define CORRELATIONMATRIX_H
#include "pairwise_matrix.h"



class CorrelationMatrix : public QAbstractTableModel, public Pairwise::Matrix
{
   Q_OBJECT
public:
   class Pair : public Matrix::Pair
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
      Pair() {}
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
   virtual QAbstractTableModel* getModel() override final { return this; }
   virtual QVariant headerData(int section, Qt::Orientation orientation, int role) const;
   virtual int rowCount(const QModelIndex&) const override final { return geneSize(); }
   virtual int columnCount(const QModelIndex&) const override final { return geneSize(); }
   virtual QVariant data(const QModelIndex& index, int role) const override final;
   void initialize(const EMetadata& geneNames, int maxClusterSize, const EMetadata& correlationNames);
   const EMetadata& correlationNames() const;
   QVector<float> dumpRawData() const;
private:
   virtual void writeHeader() { stream() << _correlationSize; }
   virtual void readHeader() { stream() >> _correlationSize; }
   static const int DATA_OFFSET {1};
   qint8 _correlationSize {0};
};



#endif
