#ifndef CORRELATIONMATRIX_H
#define CORRELATIONMATRIX_H
#include <ace/core/AceCore.h>



class CorrelationMatrix : public QAbstractTableModel, public EAbstractData
{
   Q_OBJECT
public:
   using Correlation = float;
   class Pair
   {
   public:
      Pair(CorrelationMatrix* matrix);
      Pair(const CorrelationMatrix* matrix);
      void write(int x, int y);
      void read(int x, int y) const;
      void setModeSize(qint8 size) { _modeSize = size; }
      qint8 getModeSize() const { return _modeSize; }
      Correlation& at(int mode, int index);
      const Correlation& at(int mode, int index) const;
      qint8& mode(int mode, int index);
      const qint8& mode(int mode, int index) const;
   private:
      CorrelationMatrix* _matrix {nullptr};
      const CorrelationMatrix* _cMatrix;
      mutable qint8 _modeSize;
      Correlation* _correlations;
      qint8* _modes;
   };
   virtual void readData() override final;
   virtual quint64 getDataEnd() const override final;
   virtual void newData() override final;
   virtual void prepare(bool preAllocate) override final;
   virtual void finish() override final { newData(); }
   virtual QAbstractTableModel* getModel() override final { return this; }
   virtual QVariant headerData(int section, Qt::Orientation orientation, int role) const;
   virtual int rowCount(const QModelIndex& parent) const override final;
   virtual int columnCount(const QModelIndex& parent) const override final;
   virtual QVariant data(const QModelIndex& index, int role) const override final;
   void initialize(const EMetadata& geneNames, qint32 sampleSize, qint8 correlationSize
                   , qint8 maxModes);
   qint32 getGeneSize() const { return _geneSize; }
   qint32 getSampleSize() const { return _sampleSize; }
   qint8 getCorrelationSize() const { return _correlationSize; }
   qint8 getMaxModes() const { return _maxModes; }
   static void increment(int& x, int& y);
private:
   void readPair(int x, int y, Correlation* correlations, qint8* modes, qint8& modeSize) const;
   void writePair(int x, int y, const Correlation* correlations, const qint8* modes
                  , qint8 modeSize);
   static const int DATA_OFFSET {10};
   qint32 _geneSize {0};
   qint32 _sampleSize {0};
   qint8 _correlationSize {0};
   qint8 _maxModes {0};
};



#endif
