#ifndef CORRELATIONMATRIX_H
#define CORRELATIONMATRIX_H
#include "pairwise_matrix.h"



class CorrelationMatrix : public Pairwise::Matrix
{
   Q_OBJECT
public:
   class Pair;
public:
   virtual QAbstractTableModel* model() override final;
public:
   void initialize(const EMetadata& geneNames, int maxClusterSize, const EMetadata& correlationNames);
   EMetadata correlationNames() const;
   QVector<float> dumpRawData() const;
private:
   class Model;
private:
   virtual void writeHeader() { stream() << _correlationSize; }
   virtual void readHeader() { stream() >> _correlationSize; }
   static const int DATA_OFFSET {1};
   qint8 _correlationSize {0};
  Model* _model {nullptr};
};



#endif
