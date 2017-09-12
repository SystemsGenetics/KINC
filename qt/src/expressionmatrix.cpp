#include <ace/core/metadata.h>

#include "expressionmatrix.h"






void ExpressionMatrix::readData()
{
   // seek to beginning to data
   if ( !seek(0) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed calling seek() on data object file."));
      throw e;
   }

   // read header
   stream() >> _geneSize >> _sampleSize;
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed reading from data object file."));
      throw e;
   }
}






quint64 ExpressionMatrix::getDataEnd() const
{
   // calculate and return end of data
   qint64 geneSize {_geneSize};
   qint64 sampleSize {_sampleSize};
   return DATA_OFFSET+(geneSize*sampleSize*sizeof(Expression));
}






void ExpressionMatrix::newData()
{
   // seek to beginning of data
   if ( !seek(0) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed calling seek() on data object file."));
      throw e;
   }

   // write gene and sample sizes of 0
   stream() << _geneSize << _sampleSize;
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed writing to data object file."));
      throw e;
   }
}






void ExpressionMatrix::prepare(bool preAllocate)
{
   // check if pre-allocation is requested
   if ( preAllocate )
   {
      // seek to beginning of data
      if ( !seek(0) )
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("File IO Error"));
         e.setDetails(tr("Failed calling seek() on data object file."));
         throw e;
      }

      // allocate total size needed for data
      if ( !allocate(getDataEnd()) )
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("File IO Error"));
         e.setDetails(tr("Failed allocating space in data object file."));
         throw e;
      }
   }
}






QVariant ExpressionMatrix::headerData(int section, Qt::Orientation orientation, int role) const
{
   // if this is not display role return nothing
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // get metadata root and figure out orientation
   const EMetadata::Map* map {meta().toObject()};
   switch (orientation)
   {
   case Qt::Vertical:
      // make sure metadata contains gene names
      if ( map->contains("genes") )
      {
         // get gene names and make sure it is array
         EMetadata* genes {(*map)["genes"]};
         if ( genes->isArray() )
         {
            // make sure section is within limits of array
            if ( section >= 0 && section < genes->toArray()->size() )
            {
               // return gene name
               return genes->toArray()->at(section)->toVariant();
            }
         }
      }

      // if no gene name found return nothing
      return QVariant();
   case Qt::Horizontal:
      // make sure metadata contains sample names
      if ( map->contains("samples") )
      {
         // get sample names and make sure it is array
         EMetadata* samples {(*map)["samples"]};
         if ( samples->isArray() )
         {
            // make sure section is within limits of array
            if ( section >= 0 && section < samples->toArray()->size() )
            {
               // return sample name
               return samples->toArray()->at(section)->toVariant();
            }
         }
      }

      // if no sample name found return nothing
      return QVariant();
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
   if ( !seek(DATA_OFFSET+(((index.row()*_sampleSize)+index.column())*sizeof(Expression))) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed calling seek() on data object file."));
      throw e;
   }

   // read expression from file
   stream() >> value;
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed reading from data object file."));
      throw e;
   }

   // return expression
   return value;
}






void ExpressionMatrix::initialize(QStringList geneNames, QStringList sampleNames)
{
   // get metadata root of data object
   EMetadata::Map* map {meta().toObject()};

   // create new metadata array that stores gene names and populate it
   EMetadata* metaGeneNames {new EMetadata(EMetadata::Array)};
   for (auto i = geneNames.constBegin(); i != geneNames.constEnd() ;++i)
   {
      EMetadata* name {new EMetadata(EMetadata::String)};
      *(name->toString()) = *i;
      metaGeneNames->toArray()->append(name);
   }

   // create new metadata array that stores sample names and populate it
   EMetadata* metaSampleNames {new EMetadata(EMetadata::Array)};
   for (auto i = sampleNames.constBegin(); i != sampleNames.constEnd() ;++i)
   {
      EMetadata* name {new EMetadata(EMetadata::String)};
      *(name->toString()) = *i;
      metaSampleNames->toArray()->append(name);
   }

   // insert both gene and sample names to data object's metadata
   map->insert("genes",metaGeneNames);
   map->insert("samples",metaSampleNames);

   // set gene and sample size from size of name lists
   _geneSize = geneNames.size();
   _sampleSize = sampleNames.size();
}






void ExpressionMatrix::setTransform(ExpressionMatrix::Transform transform)
{
   // get metadata root of data object and create new string
   EMetadata::Map* map {meta().toObject()};
   EMetadata* meta {new EMetadata(EMetadata::String)};

   // set new transform string according to transform type
   switch (transform)
   {
   case Transform::None:
      *(meta->toString()) = tr("none");
      break;
   case Transform::NLog:
      *(meta->toString()) = tr("natural logarithm");
      break;
   case Transform::Log2:
      *(meta->toString()) = tr("logarithm base 2");
      break;
   case Transform::Log10:
      *(meta->toString()) = tr("logarithm base 10");
      break;
   }

   // if transform key already exists in metadata remove it
   if ( map->contains("transform") )
   {
      delete map->take("transform");
   }

   // add new transform key to metadata
   map->insert("transform",meta);
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






void ExpressionMatrix::readGene(int index, Expression* expressions) const
{
   // seek to position of beginning of gene's expressions
   if ( !seek(DATA_OFFSET+(index*_sampleSize*sizeof(Expression))) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed calling seek() on data object file."));
      throw e;
   }

   // read in all expressions for gene as block of floats
   stream().read(expressions,_sampleSize);
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed reading from data object file."));
      throw e;
   }
}






void ExpressionMatrix::writeGene(int index, const Expression* expressions)
{
   // seek to position of beginning of gene's expressions
   if ( !seek(DATA_OFFSET+(index*_sampleSize*sizeof(Expression))) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed calling seek() on data object file."));
      throw e;
   }

   // overwrite all expressions for gene as block of floats
   stream().write(expressions,_sampleSize);
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed writing to data object file."));
      throw e;
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
