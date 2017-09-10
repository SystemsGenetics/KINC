#ifndef EXPRESSIONMATRIX_H
#define EXPRESSIONMATRIX_H
#include <ace/core/AceCore.h>



class ExpressionMatrix : public QAbstractTableModel, public EAbstractData
{
   Q_OBJECT
public:
   enum class Transform
   {
      None
      ,NLog
      ,Log2
      ,Log10
   };
   class Gene
   {
   public:
      Gene(ExpressionMatrix* matrix): _expressions(new float[matrix->_sampleSize]), _matrix(matrix)
         {}
      Gene(const Gene&) = delete;
      ~Gene() { delete _expressions; }
      void readGene(int index) const;
      void writeGene(int index);
      float& at(int index);
      const float& at(int index) const { return at(index); }
      float& operator[](int index) { return _expressions[index]; }
   private:
      float* _expressions;
      ExpressionMatrix* _matrix;
   };
   virtual void readData() override final;
   virtual quint64 getDataEnd() const override final;
   virtual void newData() override final;
   virtual void prepare(bool preAllocate) override final;
   virtual void finish() override final;
   virtual QAbstractTableModel* getModel() override final { return this; }
   void initialize(QStringList geneNames, QStringList sampleNames);
   void setTransform(Transform scale);
   qint32 getGeneSize() const { return _geneSize; }
   qint32 getSampleSize() const { return _sampleSize; }
   qint64 getRawSize() const;
   float* dumpRawData() const;
   QVariant headerData(int section, Qt::Orientation orientation, int role) const;
   int rowCount(const QModelIndex& parent) const override final;
   int columnCount(const QModelIndex& parent) const override final;
   QVariant data(const QModelIndex& index, int role) const override final;
private:
   void readGene(int index, float* _expressions) const;
   void writeGene(int index, const float* _expressions);
   qint32 _geneSize {0};
   qint32 _sampleSize {0};
};



#endif
