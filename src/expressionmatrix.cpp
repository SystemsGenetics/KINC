#include <ace/core/emetaarray.h>
#include <ace/core/emetaobject.h>

#include "expressionmatrix.h"






void ExpressionMatrix::readData()
{
   // seek to beginning of data
   seek(0);

   // read header
   stream() >> _geneSize >> _sampleSize;
}






qint64 ExpressionMatrix::dataEnd() const
{
   // calculate and return end of data
   qint64 geneSize {_geneSize};
   qint64 sampleSize {_sampleSize};
   return DATA_OFFSET+(geneSize*sampleSize*sizeof(Expression));
}






void ExpressionMatrix::writeNewData()
{
   // seek to beginning of data
   seek(0);

   // write gene and sample sizes of 0
   stream() << _geneSize << _sampleSize;
}






void ExpressionMatrix::finish()
{
   writeNewData();
}






QAbstractTableModel* ExpressionMatrix::model()
{
   return nullptr;
}






QVariant ExpressionMatrix::headerData(int section, Qt::Orientation orientation, int role) const
{
   // if this is not display role return nothing
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // get metadata root and figure out orientation
   const EMetaObject& object {meta().toObject()};
   switch (orientation)
   {
   case Qt::Vertical:
   {
      // get gene names and make sure it is array
      const EMetadata& genes {object["genes"]};
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
      const EMetadata& samples {object["samples"]};
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
   // get metadata root of data object
   EMetaObject& object {meta().toObject()};

   // create new metadata array that stores gene names and populate it
   EMetaArray metaGeneNames;
   for (auto i = geneNames.constBegin(); i != geneNames.constEnd() ;++i)
   {
      metaGeneNames.append(*i);
   }

   // create new metadata array that stores sample names and populate it
   EMetaArray metaSampleNames;
   for (auto i = sampleNames.constBegin(); i != sampleNames.constEnd() ;++i)
   {
      metaSampleNames.append(*i);
   }

   // insert both gene and sample names to data object's metadata
   object.insert("genes",metaGeneNames);
   object.insert("samples",metaSampleNames);

   // set gene and sample size from size of name lists
   _geneSize = geneNames.size();
   _sampleSize = sampleNames.size();
}






ExpressionMatrix::Transform ExpressionMatrix::getTransform() const
{
   // get metadata root of data object
   const EMetaObject& object {meta().toObject()};
   QString transformName = object["transform"].toString();

   // get transform type according to transform string
   if ( transformName == "none" ) {
      return Transform::None;
   }
   else if ( transformName == "natural logarithm" ) {
      return Transform::NLog;
   }
   else if ( transformName == "logarithm base 2" ) {
      return Transform::Log2;
   }
   else if ( transformName == "logarithm base 10" ) {
      return Transform::Log10;
   }

   return Transform::None;
}






void ExpressionMatrix::setTransform(ExpressionMatrix::Transform transform)
{
   // get metadata root of data object and create new string
   EMetaObject& object {meta().toObject()};
   EMetadata meta;

   // set new transform string according to transform type
   switch (transform)
   {
   case Transform::None:
      meta = tr("none");
      break;
   case Transform::NLog:
      meta = tr("natural logarithm");
      break;
   case Transform::Log2:
      meta = tr("logarithm base 2");
      break;
   case Transform::Log10:
      meta = tr("logarithm base 10");
      break;
   }

   // add new transform key to metadata
   object.insert("transform",meta);
}






qint64 ExpressionMatrix::getRawSize() const
{
   // calculate total number of floating point epxressions and return it
   qint64 geneSize {_geneSize};
   qint64 sampleSize {_sampleSize};
   return geneSize*sampleSize;
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
