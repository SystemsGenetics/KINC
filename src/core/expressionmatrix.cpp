#include "expressionmatrix.h"
#include "expressionmatrix_model.h"
//



/*!
 */
const QStringList ExpressionMatrix::_transformNames
{
   "none"
   ,"natural logarithm"
   ,"logarithm base 2"
   ,"logarithm base 10"
};
/*!
 */
const qint64 ExpressionMatrix::_dataOffset {8};






/*!
 */
qint64 ExpressionMatrix::dataEnd() const
{
   // 
   return _dataOffset + (qint64)_geneSize*(qint64)_sampleSize*sizeof(float);
}






/*!
 */
void ExpressionMatrix::readData()
{
   // 
   seek(0);
   stream() >> _geneSize >> _sampleSize;
}






/*!
 */
void ExpressionMatrix::writeNewData()
{
   // 
   setMeta(EMetadata(EMetadata::Object));

   // 
   seek(0);
   stream() << _geneSize << _sampleSize;
}






/*!
 */
void ExpressionMatrix::finish()
{
   // 
   seek(0);
   stream() << _geneSize << _sampleSize;
}






/*!
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
 */
QString ExpressionMatrix::transformString() const
{
   return meta().toObject().at("transform").toString();
}






/*!
 */
qint32 ExpressionMatrix::geneSize() const
{
   return _geneSize;
}






/*!
 */
qint32 ExpressionMatrix::sampleSize() const
{
   return _sampleSize;
}






/*!
 */
EMetadata ExpressionMatrix::geneNames() const
{
   return meta().toObject().at("genes");
}






/*!
 */
EMetadata ExpressionMatrix::sampleNames() const
{
   return meta().toObject().at("samples");
}






/*!
 */
QVector<float> ExpressionMatrix::dumpRawData() const
{
   // 
   if ( _geneSize == 0 )
   {
      return QVector<float>();
   }

   // 
   QVector<float> ret(_geneSize*_sampleSize);
   seekExpression(0,0);
   for (float& sample: ret)
   {
      stream() >> sample;
   }

   // 
   return ret;
}






/*!
 *
 * @param geneNames  
 *
 * @param sampleNames  
 *
 * @param transform  
 */
void ExpressionMatrix::initialize(const QStringList& geneNames, const QStringList& sampleNames, Transform transform)
{
   // 
   EMetaArray metaGeneNames;
   for ( auto& geneName : geneNames )
   {
      metaGeneNames.append(geneName);
   }

   // 
   EMetaArray metaSampleNames;
   for ( auto& sampleName : sampleNames )
   {
      metaSampleNames.append(sampleName);
   }

   // 
   EMetaObject metaObject {meta().toObject()};
   metaObject.insert("genes",metaGeneNames);
   metaObject.insert("samples",metaSampleNames);
   metaObject.insert("transform",_transformNames.at(static_cast<int>(transform)));
   setMeta(metaObject);

   // 
   _geneSize = geneNames.size();
   _sampleSize = sampleNames.size();
}






/*!
 *
 * @param gene  
 *
 * @param sample  
 */
void ExpressionMatrix::seekExpression(int gene, int sample) const
{
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
   seek(_dataOffset + ((qint64)gene*(qint64)_sampleSize + (qint64)sample)*sizeof(float));
}
