#ifndef CCMATRIX_H
#define CCMATRIX_H
#include "pairwise_matrix.h"



class CCMatrix : public Pairwise::Matrix
{
   Q_OBJECT
public:
   class Pair;
public:
   virtual QAbstractTableModel* model() override final;
public:
   void initialize(const EMetadata& geneNames, int maxClusterSize, const EMetadata& sampleNames);
   EMetadata sampleNames() const;
   int sampleSize() const { return _sampleSize; }
private:
   class Model;
private:
   virtual void writeHeader() { stream() << _sampleSize; }
   virtual void readHeader() { stream() >> _sampleSize; }
   static const int DATA_OFFSET {4};
   qint32 _sampleSize {0};
  Model* _model {nullptr};
};



#endif
