#ifndef CORRELATIONMATRIX_H
#define CORRELATIONMATRIX_H
#include "pairwise_matrix.h"



class CorrelationMatrix : public Pairwise::Matrix
{
   Q_OBJECT
public:
   class Pair;
   virtual QAbstractTableModel* model() override final;
   QVariant headerData(int section, Qt::Orientation orientation, int role) const;
   int rowCount(const QModelIndex&) const;
   int columnCount(const QModelIndex&) const;
   QVariant data(const QModelIndex& index, int role) const;
   void initialize(const EMetadata& geneNames, int maxClusterSize, const EMetadata& correlationNames);
   EMetadata correlationNames() const;
   QVector<float> dumpRawData() const;
private:
   virtual void writeHeader() { stream() << _correlationSize; }
   virtual void readHeader() { stream() >> _correlationSize; }
   static const int DATA_OFFSET {1};
   qint8 _correlationSize {0};
};



#endif
