#include "expressionmatrix.h"

#define DATA_OFFSET 8






void ExpressionMatrix::readData()
{
}






quint64 ExpressionMatrix::getDataEnd() const
{
   qint64 geneSize {_geneSize};
   qint64 sampleSize {_sampleSize};
   return DATA_OFFSET+(geneSize*sampleSize);
}






void ExpressionMatrix::newData()
{
}






void ExpressionMatrix::prepare(bool preAllocate)
{
}






void ExpressionMatrix::finish()
{
}






void ExpressionMatrix::initialize(QStringList geneNames, QStringList sampleNames)
{
}






QVariant ExpressionMatrix::headerData(int section, Qt::Orientation orientation, int role) const
{
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
