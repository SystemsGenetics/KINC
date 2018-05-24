#include "expressionmatrix.h"






const QStringList ExpressionMatrix::TRANSFORM_NAMES
{
   "none"
   ,"natural logarithm"
   ,"logarithm base 2"
   ,"logarithm base 10"
};






qint64 ExpressionMatrix::dataEnd() const
{
   // calculate and return end of data
   return DATA_OFFSET + ((qint64)_geneSize * (qint64)_sampleSize * sizeof(Expression));
}






void ExpressionMatrix::readData()
{
   // read header
   seek(0);
   stream() >> _geneSize >> _sampleSize;
}






void ExpressionMatrix::writeNewData()
{
   // initialize metadata
   setMeta(EMetadata(EMetadata::Object));

   // initialize header
   seek(0);
   stream() << _geneSize << _sampleSize;
}






QAbstractTableModel* ExpressionMatrix::model()
{
   return nullptr;
}






void ExpressionMatrix::finish()
{
   // write header
   seek(0);
   stream() << _geneSize << _sampleSize;
}






QVariant ExpressionMatrix::headerData(int section, Qt::Orientation orientation, int role) const
{
   // if this is not display role return nothing
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // get metadata root and figure out orientation
   switch (orientation)
   {
   case Qt::Vertical:
   {
      // get gene names and make sure it is array
      const EMetadata& genes {meta().toObject().at("genes")};
      if ( genes.isArray() )
      {
         // make sure section is within limits of array
         if ( section >= 0 && section < genes.toArray().size() )
         {
            // return gene name
            return genes.toArray().at(section).toString();
         }
      }

      // if no gene name found return nothing
      return QVariant();
   }
   case Qt::Horizontal:
   {
      // get sample names and make sure it is array
      const EMetadata& samples {meta().toObject().at("samples")};
      if ( samples.isArray() )
      {
         // make sure section is within limits of array
         if ( section >= 0 && section < samples.toArray().size() )
         {
            // return sample name
            return samples.toArray().at(section).toString();
         }
      }

      // if no sample name found return nothing
      return QVariant();
   }
   default:
      // unknown orientation so return nothing
      return QVariant();
   }
}






int ExpressionMatrix::rowCount(const QModelIndex& parent) const
{
   // return gene size for row count
   Q_UNUSED(parent);
   return _geneSize;
}






int ExpressionMatrix::columnCount(const QModelIndex& parent) const
{
   // return sample size for column count
   Q_UNUSED(parent);
   return _sampleSize;
}






QVariant ExpressionMatrix::data(const QModelIndex& index, int role) const
{
   // if role is not display return nothing
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // if index is out of range return nothing
   if ( index.row() >= _geneSize || index.column() >= _sampleSize )
   {
      return QVariant();
   }

   // make input variable and seek to position of queried expression
   Expression value;
   seek(DATA_OFFSET+(((index.row()*_sampleSize)+index.column())*sizeof(Expression)));

   // read expression from file
   stream() >> value;

   // return expression
   return value;
}






void ExpressionMatrix::initialize(QStringList geneNames, QStringList sampleNames)
{
   // create metadata array of gene names
   EMetaArray metaGeneNames;
   for ( auto& geneName : geneNames )
   {
      metaGeneNames.append(geneName);
   }

   // create metadata array of sample names
   EMetaArray metaSampleNames;
   for ( auto& sampleName : sampleNames )
   {
      metaSampleNames.append(sampleName);
   }

   // insert gene and sample names to data object's metadata
   meta().toObject().insert("genes", metaGeneNames);
   meta().toObject().insert("samples", metaSampleNames);

   // set gene and sample size
   _geneSize = geneNames.size();
   _sampleSize = sampleNames.size();
}






ExpressionMatrix::Transform ExpressionMatrix::getTransform() const
{
   auto& transformName {meta().toObject().at("transform").toString()};
   return static_cast<Transform>(TRANSFORM_NAMES.indexOf(transformName));
}






void ExpressionMatrix::setTransform(ExpressionMatrix::Transform transform)
{
   auto& transformName {TRANSFORM_NAMES.at(static_cast<int>(transform))};
   meta().toObject().insert("transform", transformName);
}






qint64 ExpressionMatrix::getRawSize() const
{
   return (qint64)_geneSize * (qint64)_sampleSize;
}






ExpressionMatrix::Expression* ExpressionMatrix::dumpRawData() const
{
   // if there are no genes do nothing
   if ( _geneSize == 0 )
   {
      return nullptr;
   }

   // create new floating point array and populate with all gene expressions
   Expression* ret {new Expression[getRawSize()]};
   for (int i = 0; i < _geneSize ;++i)
   {
      readGene(i,&ret[i*_sampleSize]);
   }

   // return new float array
   return ret;
}






EMetadata ExpressionMatrix::getGeneNames() const
{
   return meta().toObject().at("genes");
}






EMetadata ExpressionMatrix::getSampleNames() const
{
   return meta().toObject().at("samples");
}






void ExpressionMatrix::readGene(int index, Expression* expressions) const
{
   // seek to position of beginning of gene's expressions
   seek(DATA_OFFSET + (index * _sampleSize * sizeof(Expression)));

   // read in all expressions for gene as block of floats
   for ( int i = 0; i < _sampleSize; ++i )
   {
      stream() >> expressions[i];
   }
}






void ExpressionMatrix::writeGene(int index, const Expression* expressions)
{
   // seek to position of beginning of gene's expressions
   seek(DATA_OFFSET + (index * _sampleSize * sizeof(Expression)));

   // overwrite all expressions for gene as block of floats
   for ( int i = 0; i < _sampleSize; ++i )
   {
      stream() << expressions[i];
   }
}






void ExpressionMatrix::Gene::read(int index) const
{
   // make sure given gene index is within range
   if ( index < 0 || index >= _matrix->_geneSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to read gene %1 when maximum is %2.").arg(index)
                   .arg(_matrix->_geneSize-1));
      throw e;
   }

   // read gene expressions from data object
   _matrix->readGene(index,_expressions);
}






void ExpressionMatrix::Gene::write(int index)
{
   // make sure given gene index is within range
   if ( index < 0 || index >= _matrix->_geneSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to write gene %1 when maximum is %2.").arg(index)
                   .arg(_matrix->_geneSize-1));
      throw e;
   }

   // write gene expressions to data object
   _matrix->writeGene(index,_expressions);
}






ExpressionMatrix::Expression& ExpressionMatrix::Gene::at(int index)
{
   // make sure given sample index is within range
   if ( index < 0 || index >= _matrix->_sampleSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to access gene expression %1 when maximum is %2.").arg(index)
                   .arg(_matrix->_sampleSize-1));
      throw e;
   }

   // return gene expression
   return _expressions[index];
}
