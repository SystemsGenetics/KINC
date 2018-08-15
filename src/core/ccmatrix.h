#ifndef CCMATRIX_H
#define CCMATRIX_H
#include "pairwise_matrix.h"



class CCMatrix : public Pairwise::Matrix
{
   Q_OBJECT
public:
   class Pair;
   virtual QAbstractTableModel* model() override final;
   QVariant headerData(int section, Qt::Orientation orientation, int role) const;
   int rowCount(const QModelIndex&) const;
   int columnCount(const QModelIndex&) const;
   QVariant data(const QModelIndex& index, int role) const;
   void initialize(const EMetadata& geneNames, int maxClusterSize, const EMetadata& sampleNames);
   EMetadata sampleNames() const;
   int sampleSize() const { return _sampleSize; }
private:
   virtual void writeHeader() { stream() << _sampleSize; }
   virtual void readHeader() { stream() >> _sampleSize; }
   static const int DATA_OFFSET {4};
   qint32 _sampleSize {0};
};



#endif
