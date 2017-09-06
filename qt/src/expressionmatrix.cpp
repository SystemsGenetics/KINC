#include <ace/core/metadata.h>

#include "expressionmatrix.h"

#define DATA_OFFSET 8






void ExpressionMatrix::readData()
{
   if ( !seek(0) )
   {
      ;//ERROR
   }
   stream() >> _geneSize >> _sampleSize;
   if ( !stream() )
   {
      ;//ERRPOR
   }
}






quint64 ExpressionMatrix::getDataEnd() const
{
   qint64 geneSize {_geneSize};
   qint64 sampleSize {_sampleSize};
   return DATA_OFFSET+(geneSize*sampleSize);
}






void ExpressionMatrix::newData()
{
   if ( !seek(0) )
   {
      ;//ERROR
   }
   stream() << _geneSize << _sampleSize;
   if ( !stream() )
   {
      ;//ERRPOR
   }
}






void ExpressionMatrix::prepare(bool preAllocate)
{
   if ( preAllocate )
   {
      if ( !seek(0) )
      {
         ;//ERROR
      }
      if ( !allocate(DATA_OFFSET+(_geneSize*_sampleSize*sizeof(float))) )
      {
         ;//ERROR
      }
   }
}






void ExpressionMatrix::finish()
{
   if ( !seek(0) )
   {
      ;//ERROR
   }
   stream() << _geneSize << _sampleSize;
   if ( !stream() )
   {
      ;//ERRPOR
   }
}






void ExpressionMatrix::initialize(QStringList geneNames, QStringList sampleNames)
{
   EMetadata::Map* map {meta().toObject()};
   EMetadata* metaGeneNames {new EMetadata(EMetadata::Array)};
   for (auto i = geneNames.constBegin(); i != geneNames.constEnd() ;++i)
   {
      EMetadata* name {new EMetadata(EMetadata::String)};
      *(name->toString()) = *i;
      metaGeneNames->toArray()->append(name);
   }
   EMetadata* metaSampleNames {new EMetadata(EMetadata::Array)};
   for (auto i = sampleNames.constBegin(); i != sampleNames.constEnd() ;++i)
   {
      EMetadata* name {new EMetadata(EMetadata::String)};
      *(name->toString()) = *i;
      metaSampleNames->toArray()->append(name);
   }
   map->insert("genes",metaGeneNames);
   map->insert("samples",metaSampleNames);
   _geneSize = geneNames.size();
   _sampleSize = sampleNames.size();
}






QVariant ExpressionMatrix::headerData(int section, Qt::Orientation orientation, int role) const
{
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }
   EMetadata::Map* map {meta().toObject()};
   switch (orientation)
   {
   case Qt::Vertical:
      if ( map->contains("genes") )
      {
         EMetadata* genes {(*map)["genes"]};
         if ( genes->isArray() )
         {
            if ( section >= 0 && section < genes->toArray()->size() )
            {
               return genes->toArray()->at(section)->toVariant();
            }
         }
      }
   case Qt::Horizontal:
      break;
   default:
      return QVariant();
   }
}






int ExpressionMatrix::rowCount(const QModelIndex& parent) const
{
}






int ExpressionMatrix::columnCount(const QModelIndex& parent) const
{
}






QVariant ExpressionMatrix::data(const QModelIndex& index, int role) const
{
}






void ExpressionMatrix::syncGene(int index, float* expressions, Sync which)
{
   if ( _newData )
   {
      ;//ERROR
   }
   if ( !seek(DATA_OFFSET+(_sampleSize*index)) )
   {
      ;//ERROR
   }
   for (int i = 0; i < _sampleSize ;++i)
   {
      switch (which)
      {
      case Sync::Read:
         stream() >> expressions[i];
         break;
      case Sync::Write:
         stream() << expressions[i];
         break;
      }
   }
   if ( !stream() )
   {
      ;//ERROR!
   }
}






void ExpressionMatrix::Gene::readGene(int index)
{
   if ( index < 0 || index >= _matrix->_geneSize)
   {
      ;//ERROR
   }
   _matrix->syncGene(index,_expressions,Sync::Read);
}






void ExpressionMatrix::Gene::writeGene(int index) const
{
   if ( index < 0 || index >= _matrix->_geneSize)
   {
      ;//ERROR
   }
   _matrix->syncGene(index,_expressions,Sync::Write);
}






float& ExpressionMatrix::Gene::at(int index)
{
   if ( index < 0 || index >= _matrix->_sampleSize )
   {
      ;//ERROR
   }
   return _expressions[index];
}
