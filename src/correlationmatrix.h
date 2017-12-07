#ifndef CORRELATIONMATRIX_H
#define CORRELATIONMATRIX_H
#include "genepair_base.h"



class CorrelationMatrix : public QAbstractTableModel, public GenePair::Base
{
   Q_OBJECT
public:
   using GenePair::Vector;
   class Pair
   {
   public:
      Pair(CorrelationMatrix* matrix);
      Pair(const CorrelationMatrix* matrix);
      void write(Vector index) { _matrix->writePair(index,_correlations); }
      void read(Vector index) const { _correlations = _cMatrix->readPair(index); }
      void readFirst() const;
      bool readNext() const;
      void clear() const;
      int addCluster() const;
      int clusterSize() const;
      bool isEmpty() const;
      const float& at(int cluster, int correlation) const;
   private:
      CorrelationMatrix* _matrix {nullptr};
      const CorrelationMatrix* _cMatrix;
      mutable QList<QList<float>> _correlations;
      mutable qint64 _index {-1};
   };
   virtual void readData() override final;
   virtual void finish() override final;
   virtual QAbstractTableModel* getModel() override final { return this; }
   virtual QVariant headerData(int section, Qt::Orientation orientation, int role) const;
   virtual int rowCount(const QModelIndex&) const override final { return _geneSize; }
   virtual int columnCount(const QModelIndex&) const override final { return _geneSize; }
   virtual QVariant data(const QModelIndex& index, int role) const override final;
   void initialize(const EMetadata& geneNames, qint32 sampleSize
                   , const EMetadata& correlationNames);
   qint32 geneSize() const { return _geneSize; }
   qint32 sampleSize() const { return _sampleSize; }
   qint8 correlationSize() const { return _correlationSize; }
   const EMetadata& getGeneNames() const;
private:
   QList<QList<float>> readPair(Vector index) const;
   QList<QList<float>> readPair(qint64* index) const;
   void writePair(Vector index, const QList<QList<float>>& correlations);
   static const int DATA_OFFSET {9};
   qint32 _geneSize {0};
   qint32 _sampleSize {0};
   qint8 _correlationSize {0};
};



#endif
