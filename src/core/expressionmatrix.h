#ifndef EXPRESSIONMATRIX_H
#define EXPRESSIONMATRIX_H
#include <ace/core/core.h>
//



/*!
 * This class implements the expression matrix data object. An expression matrix
 * is a matrix of real numbers whose rows represent genes and whose columns
 * represent samples. The matrix data can be accessed using the gene interator,
 * which iterates through each gene (row) in the matrix.
 */
class ExpressionMatrix : public EAbstractData
{
   Q_OBJECT
public:
   class Gene;
public:
   virtual qint64 dataEnd() const override final;
   virtual void readData() override final;
   virtual void writeNewData() override final;
   virtual void finish() override final;
   virtual QAbstractTableModel* model() override final;
public:
   qint32 geneSize() const;
   qint32 sampleSize() const;
   EMetaArray geneNames() const;
   EMetaArray sampleNames() const;
   std::vector<float> dumpRawData() const;
   void initialize(const QStringList& geneNames, const QStringList& sampleNames);
private:
   class Model;
private:
   void seekExpression(int gene, int sample) const;
   /*!
    * The header size (in bytes) at the beginning of the file. The header
    * consists of the gene size and the sample size.
    */
   constexpr static const qint64 _headerSize {8};
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
