#ifndef EXPRESSIONMATRIX_H
#define EXPRESSIONMATRIX_H
#include <ace/core/core.h>



class ExpressionMatrix : public EAbstractData
{
   Q_OBJECT
public:
   using Expression = float;
   static const QStringList TRANSFORM_NAMES;
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
      Gene(ExpressionMatrix* matrix):
         _expressions(new Expression[matrix->_sampleSize]),
         _matrix(matrix)
         {}
      Gene(const Gene&) = delete;
      ~Gene() { delete _expressions; }
      void read(int index) const;
      void write(int index);
      Expression& at(int index);
      const Expression& at(int index) const;
      Expression& operator[](int index) { return _expressions[index]; }
   private:
      Expression* _expressions;
      ExpressionMatrix* _matrix;
   };
   virtual qint64 dataEnd() const override final;
   virtual void readData() override final;
   virtual void writeNewData() override final;
   virtual QAbstractTableModel* model() override final;
   virtual void finish() override final;
   QVariant headerData(int section, Qt::Orientation orientation, int role) const;
   int rowCount(const QModelIndex& parent) const;
   int columnCount(const QModelIndex& parent) const;
   QVariant data(const QModelIndex& index, int role) const;
   void initialize(QStringList geneNames, QStringList sampleNames);
   Transform getTransform() const;
   void setTransform(Transform scale);
   qint32 getGeneSize() const { return _geneSize; }
   qint32 getSampleSize() const { return _sampleSize; }
   qint64 getRawSize() const;
   Expression* dumpRawData() const;
   EMetadata getGeneNames() const;
   EMetadata getSampleNames() const;
private:
   void readGene(int index, Expression* expressions) const;
   void writeGene(int index, const Expression* expressions);
   static const int DATA_OFFSET {8};
   qint32 _geneSize {0};
   qint32 _sampleSize {0};
};



#endif
