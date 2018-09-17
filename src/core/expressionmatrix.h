#ifndef EXPRESSIONMATRIX_H
#define EXPRESSIONMATRIX_H
#include <ace/core/core.h>
//



/*!
 * This class implements the expression matrix data object. An expression matrix
 * is a matrix of real numbers whose rows represent genes and whose columns
 * represent samples. The matrix data can be accessed using the Gene interator,
 * which iterates through each gene (row) in the matrix.
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
    * Defines the transforms that this data object supports.
    */
   enum class Transform
   {
      /*!
       * no transform
       */
      None
      /*!
       * natural logarithm
       */
      ,NatLog
      /*!
       * base-2 logarithm
       */
      ,Log2
      /*!
       * base-10 logarithm
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
    * String list of transforms for this data object corresponding to the
    * Transform type.
    */
   static const QStringList _transformNames;
private:
   class Model;
private:
   void seekExpression(int gene, int sample) const;
   /*!
    * The header size (in bytes) at the beginning of the file.
    */
   static const qint64 _dataOffset;
   /*!
    * The number of genes (rows) in the expression matrix.
    */
   qint32 _geneSize;
   /*!
    * The number of samples (columns) in the expression matrix.
    */
   qint32 _sampleSize;
   /*!
    * Pointer to a qt table model for this class.
    */
   Model* _model {nullptr};
};



#endif
