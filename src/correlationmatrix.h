#ifndef CORRELATIONMATRIX_H
#define CORRELATIONMATRIX_H
#include "correlationbase.h"



class CorrelationMatrix : public QAbstractTableModel, public CorrelationBase
{
   Q_OBJECT
public:
   class Pair
   {
   public:
      Pair(CorrelationMatrix* matrix);
      Pair(const CorrelationMatrix* matrix);
      void write(Iterator index) { _matrix->writePair(index,_correlations); }
      void read(Iterator index) const { _cMatrix->readPair(index,_correlations); }
      void readFirst() const;
      bool readNext() const;
      float& at(int mode, int correlation);
      const float& at(int mode, int correlation) const;
   private:
      CorrelationMatrix* _matrix {nullptr};
      const CorrelationMatrix* _cMatrix;
      float* _correlations;
      mutable qint64 _index {-1};
   };
   virtual void readData() override final;
   virtual void finish() override final;
   virtual QAbstractTableModel* getModel() override final { return this; }
   virtual QVariant headerData(int section, Qt::Orientation orientation, int role) const;
   virtual int rowCount(const QModelIndex& parent) const override final;
   virtual int columnCount(const QModelIndex& parent) const override final;
   virtual QVariant data(const QModelIndex& index, int role) const override final;
   void initialize(const EMetadata& geneNames, qint32 sampleSize
                   , const EMetadata& correlationNames, qint8 maxModes);
   qint32 getGeneSize() const { return _geneSize; }
   qint32 getSampleSize() const { return _sampleSize; }
   qint8 getCorrelationSize() const { return _correlationSize; }
   qint8 getMaxModes() const { return _maxModes; }
   const EMetadata& getGeneNames() const;
private:
   void readPair(Iterator index, float* correlations) const;
   void readPair(qint64 index, float* correlations) const;
   void writePair(Iterator index, const float* correlations);
   static const int DATA_OFFSET {10};
   qint32 _geneSize {0};
   qint32 _sampleSize {0};
   qint8 _correlationSize {0};
   qint8 _maxModes {0};
};



#endif
