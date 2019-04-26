#include "expressionmatrix.h"
#include "expressionmatrix_model.h"
//



/*!
 * Return the index of the first byte in this data object after the end of
 * the data section. Defined as the header size plus the size of the matrix data.
 */
qint64 ExpressionMatrix::dataEnd() const
{
   EDEBUG_FUNC(this);

   return _headerSize + static_cast<qint64>(_geneSize) * _sampleSize * sizeof(float);
}






/*!
 * Read in the data of an existing data object that was just opened.
 */
void ExpressionMatrix::readData()
{
   EDEBUG_FUNC(this);

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
   EDEBUG_FUNC(this);

   // initialize metadata object
   setMeta(EMetaObject());

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
   EDEBUG_FUNC(this);

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
   EDEBUG_FUNC(this);

   if ( !_model )
   {
      _model = new Model(this);
   }
   return _model;
}






/*!
 * Return the number of genes (rows) in this expression matrix.
 */
qint32 ExpressionMatrix::geneSize() const
{
   EDEBUG_FUNC(this);

   return _geneSize;
}






/*!
 * Return the number of samples (columns) in this expression matrix.
 */
qint32 ExpressionMatrix::sampleSize() const
{
   EDEBUG_FUNC(this);

   return _sampleSize;
}






/*!
 * Return the list of gene names in this expression matrix.
 */
EMetaArray ExpressionMatrix::geneNames() const
{
   EDEBUG_FUNC(this);

   return meta().toObject().at("genes").toArray();
}






/*!
 * Return the list of sample names in this expression matrix.
 */
EMetaArray ExpressionMatrix::sampleNames() const
{
   EDEBUG_FUNC(this);

   return meta().toObject().at("samples").toArray();
}






/*!
 * Return an array of this expression matrix's data in row-major order.
 */
std::vector<float> ExpressionMatrix::dumpRawData() const
{
   EDEBUG_FUNC(this);

   // return empty array if expression matrix is empty
   if ( _geneSize == 0 )
   {
      return std::vector<float>();
   }

   // allocate an array with the same size as the expression matrix
   std::vector<float> ret(static_cast<qint64>(_geneSize) * _sampleSize);

   // seek to the beginning of the expression data
   seekExpression(0,0);

   // write each expression to the array
   for ( float& sample: ret )
   {
      stream() >> sample;
   }

   // return the array
   return ret;
}






/*!
 * Initialize this expression matrix with a list of gene names and a list of
 * sample names.
 *
 * @param geneNames
 * @param sampleNames
 */
void ExpressionMatrix::initialize(const QStringList& geneNames, const QStringList& sampleNames)
{
   EDEBUG_FUNC(this,&geneNames,&sampleNames);

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

   // save the gene names and sample names to metadata
   EMetaObject metaObject {meta().toObject()};
   metaObject.insert("genes",metaGeneNames);
   metaObject.insert("samples",metaSampleNames);
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
   EDEBUG_FUNC(this,gene,sample);

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
   seek(_headerSize + (static_cast<qint64>(gene) * _sampleSize + sample) * sizeof(float));
}
