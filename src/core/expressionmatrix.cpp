#include "expressionmatrix.h"
#include "expressionmatrix_model.h"
//



/*!
 * String list of transforms for this data object corresponding to the
 * Transform type.
 */
const QStringList ExpressionMatrix::_transformNames
{
   "none"
   ,"natural logarithm"
   ,"logarithm base 2"
   ,"logarithm base 10"
};
/*!
 * The header size (in bytes) at the beginning of the file. The header
 * consists of the gene size and the sample size.
 */
const qint64 ExpressionMatrix::_dataOffset {8};






/*!
 * Return the index of the first byte in this data object after the end of
 * the data section.
 */
qint64 ExpressionMatrix::dataEnd() const
{
   return _dataOffset + (qint64)_geneSize*(qint64)_sampleSize*sizeof(float);
}






/*!
 * Read in the data of an existing data object that was just opened.
 */
void ExpressionMatrix::readData()
{
   // seek to the beginning of the data
   seek(0);

   // read the header
   stream() >> _geneSize >> _sampleSize;
}






/*!
 * Initialize this data object's data to a null state.
 */
void ExpressionMatrix::writeNewData()
{
   // initialize metadata object
   setMeta(EMetadata(EMetadata::Object));

   // seek to the beginning of the data
   seek(0);

   // write the header
   stream() << _geneSize << _sampleSize;
}






/*!
 * Finalize this data object's data after the analytic that created it has
 * finished giving it new data.
 */
void ExpressionMatrix::finish()
{
   // seek to the beginning of the data
   seek(0);

   // write the header
   stream() << _geneSize << _sampleSize;
}






/*!
 * Return a qt table model that represents this data object as a table.
 */
QAbstractTableModel* ExpressionMatrix::model()
{
   if ( !_model )
   {
      _model = new Model(this);
   }
   return _model;
}






/*!
 * Return the transform that was applied to this expression matrix.
 */
ExpressionMatrix::Transform ExpressionMatrix::transform() const
{
   QString transformName {meta().toObject().at("transform").toString()};
   int index {_transformNames.indexOf(transformName)};
   if ( index == -1 )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Logic Error"));
      e.setDetails(tr("Unknown transform name: %1.").arg(transformName));
      throw e;
   }
   return static_cast<Transform>(index);
}






/*!
 * Return the name of transform that was applied to this expression matrix.
 */
QString ExpressionMatrix::transformString() const
{
   return meta().toObject().at("transform").toString();
}






/*!
 * Return the number of genes (rows) in this expression matrix.
 */
qint32 ExpressionMatrix::geneSize() const
{
   return _geneSize;
}






/*!
 * Return the number of samples (columns) in this expression matrix.
 */
qint32 ExpressionMatrix::sampleSize() const
{
   return _sampleSize;
}






/*!
 * Return the list of gene names in this expression matrix.
 */
EMetadata ExpressionMatrix::geneNames() const
{
   return meta().toObject().at("genes");
}






/*!
 * Return the list of sample names in this expression matrix.
 */
EMetadata ExpressionMatrix::sampleNames() const
{
   return meta().toObject().at("samples");
}






/*!
 * Return an array of this expression matrix's data in row-major form.
 */
QVector<float> ExpressionMatrix::dumpRawData() const
{
   // return empty array if expression matrix is empty
   if ( _geneSize == 0 )
   {
      return QVector<float>();
   }

   // allocate an array with the same size as the expression matrix
   QVector<float> ret(_geneSize*_sampleSize);

   // seek to the beginning of the expression data
   seekExpression(0,0);

   // write each expression to the array
   for (float& sample: ret)
   {
      stream() >> sample;
   }

   // return the array
   return ret;
}






/*!
 * Initialize this expression matrix with a list of gene names, a list of
 * sample names, and a transform.
 *
 * @param geneNames
 * @param sampleNames
 * @param transform
 */
void ExpressionMatrix::initialize(const QStringList& geneNames, const QStringList& sampleNames, Transform transform)
{
   // create a metadata array of gene names
   EMetaArray metaGeneNames;
   for ( auto& geneName : geneNames )
   {
      metaGeneNames.append(geneName);
   }

   // create a metadata array of sample names
   EMetaArray metaSampleNames;
   for ( auto& sampleName : sampleNames )
   {
      metaSampleNames.append(sampleName);
   }

   // save the gene names, sample names, and transform to metadata
   EMetaObject metaObject {meta().toObject()};
   metaObject.insert("genes",metaGeneNames);
   metaObject.insert("samples",metaSampleNames);
   metaObject.insert("transform",_transformNames.at(static_cast<int>(transform)));
   setMeta(metaObject);

   // initialize the gene size and sample size accordingly
   _geneSize = geneNames.size();
   _sampleSize = sampleNames.size();
}






/*!
 * Seek to a particular expression in this expression matrix given a gene index
 * and a sample index.
 *
 * @param gene
 * @param sample
 */
void ExpressionMatrix::seekExpression(int gene, int sample) const
{
   // make sure that the indices are valid
   if ( gene < 0 || gene >= _geneSize || sample < 0 || sample >= _sampleSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Invalid Argument"));
      e.setDetails(tr("Invalid (gene,sample) index (%1,%2) with size of (%1,%2).")
                   .arg(gene)
                   .arg(sample)
                   .arg(_geneSize)
                   .arg(_sampleSize));
      throw e;
   }

   // seek to the specified position in the data
   seek(_dataOffset + ((qint64)gene*(qint64)_sampleSize + (qint64)sample)*sizeof(float));
}
