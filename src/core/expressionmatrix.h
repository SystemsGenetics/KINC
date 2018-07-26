#ifndef EXPRESSIONMATRIX_H
#define EXPRESSIONMATRIX_H
#include <ace/core/core.h>
//



/*!
 */
class ExpressionMatrix : public EAbstractData
{
   Q_OBJECT
public:
   class Gene;
   class ConstGene;
public:
   virtual qint64 dataEnd() const override final;
   virtual void readData() override final;
   virtual void writeNewData() override final;
   virtual void finish() override final;
   virtual QAbstractTableModel* model() override final;
public:
   /*!
    */
   enum class Transform
   {
      /*!
       */
      None
      /*!
       */
      ,NatLog
      /*!
       */
      ,Log2
      /*!
       */
      ,Log10
   };
   ExpressionMatrix::Transform transform() const;
   QString transformString() const;
   qint32 geneSize() const;
   qint32 sampleSize() const;
   EMetadata geneNames() const;
   EMetadata sampleNames() const;
   QVector<float> dumpRawData() const;
   void initialize(const QStringList& geneNames, const QStringList& sampleNames, Transform transform);
   /*!
    */
   static const QStringList _transformNames;
private:
   class Model;
private:
   void seekExpression(int gene, int sample) const;
   /*!
    */
   static const qint64 _dataOffset;
   /*!
    */
   qint32 _geneSize;
   /*!
    */
   qint32 _sampleSize;
   /*!
    */
   Model* _model {nullptr};
};



#endif
