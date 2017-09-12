#include <ace/core/metadata.h>

#include "correlationmatrix.h"






void CorrelationMatrix::readData()
{
   if ( !seek(0) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed calling seek() on data object file."));
      throw e;
   }
   stream() >> _geneSize >> _sampleSize >> _correlationSize >> _maxModes;
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed reading from data object file."));
      throw e;
   }
}






quint64 CorrelationMatrix::getDataEnd() const
{
   qint64 geneSize {_geneSize};
   qint64 correlationSize = _correlationSize*_maxModes*sizeof(Correlation);
   qint64 modeSize = _sampleSize*_maxModes*sizeof(qint8);
   return (geneSize*(geneSize-1)/2)*(correlationSize+modeSize+sizeof(qint8));
}






void CorrelationMatrix::newData()
{
   if ( !seek(0) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed calling seek() on data object file."));
      throw e;
   }
   stream() << _geneSize << _sampleSize << _correlationSize << _maxModes;
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed writing to data object file."));
      throw e;
   }
}






void CorrelationMatrix::prepare(bool preAllocate)
{
   if ( preAllocate )
   {
      if ( !seek(0) )
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("File IO Error"));
         e.setDetails(tr("Failed calling seek() on data object file."));
         throw e;
      }
      if ( !allocate(getDataEnd()) )
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("File IO Error"));
         e.setDetails(tr("Failed allocating space in data object file."));
         throw e;
      }
   }
}






QVariant CorrelationMatrix::headerData(int section, Qt::Orientation orientation, int role) const
{
   Q_UNUSED(orientation);
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }
   const EMetadata::Map* map {meta().toObject()};
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
   return QVariant();
}






int CorrelationMatrix::rowCount(const QModelIndex& parent) const
{
   Q_UNUSED(parent);
   return _geneSize;
}






int CorrelationMatrix::columnCount(const QModelIndex& parent) const
{
   Q_UNUSED(parent);
   return _geneSize;
}






QVariant CorrelationMatrix::data(const QModelIndex& index, int role) const
{
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }
   if ( index.row() == index.column() )
   {
      return 1;
   }
   const Pair pair(this);
   pair.read(index.row(),index.column());
   QString ret;
   for (int i = 0; i < pair.getModeSize() ;++i)
   {
      ret.append(tr("Mode %1: ").arg(i));
      for (int j = 0; j < _correlationSize ;)
      {
         ret.append(QString::number(pair.at(i,j)));
         if ( ++j < _correlationSize )
         {
            ret.append(", ");
         }
      }
      ret.append("\n");
   }
   return ret;
}






void CorrelationMatrix::initialize(EMetadata* geneNames, qint32 sampleSize, qint8 correlationSize
                                   , qint8 maxModes)
{
   if ( !geneNames->isArray() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Gene names metadata is not correct array type."));
      throw e;
   }
   EMetadata::Map* map {meta().toObject()};
   map->insert("genes",new EMetadata(*geneNames));
   _geneSize = geneNames->toArray()->size();
   _sampleSize = sampleSize;
   _correlationSize = correlationSize;
   _maxModes = maxModes;
}






void CorrelationMatrix::increment(int &x, int &y)
{
   if ( x >= ++y )
   {
      ++x;
      y = 0;
   }
}






void CorrelationMatrix::readCorrelation(int x, int y, Correlation* correlations, qint8* modes
                                        , qint8& modeSize) const
{
   qint64 bX {x};
   qint64 bY {y};
   qint64 index {(bX*(bX-1)/2)+bY};
   qint64 correlationSize = _correlationSize*_maxModes*sizeof(Correlation);
   qint64 iModeSize = _sampleSize*_maxModes*sizeof(qint8);
   if ( !seek(DATA_OFFSET+(index*(correlationSize+iModeSize+sizeof(qint8)))) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed calling seek() on data object file."));
      throw e;
   }
   stream() >> modeSize;
   stream().read(correlations,_correlationSize);
   stream().read(modes,iModeSize);
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed reading from data object file."));
      throw e;
   }
}






void CorrelationMatrix::writeCorrelation(int x, int y, const Correlation* correlations
                                         , const qint8* modes, qint8 modeSize)
{
   qint64 bX {x};
   qint64 bY {y};
   qint64 index {(bX*(bX-1)/2)+bY};
   qint64 correlationSize = _correlationSize*_maxModes*sizeof(Correlation);
   qint64 iModeSize = _sampleSize*_maxModes*sizeof(qint8);
   if ( !seek(DATA_OFFSET+(index*(correlationSize+iModeSize+sizeof(qint8)))) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed calling seek() on data object file."));
      throw e;
   }
   stream() << modeSize;
   stream().write(correlations,_correlationSize);
   stream().write(modes,iModeSize);
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed reading from data object file."));
      throw e;
   }
}






CorrelationMatrix::Pair::Pair(CorrelationMatrix* matrix):
   _matrix(matrix),
   _cMatrix(matrix)
{
   _correlations = new Correlation[_matrix->_correlationSize*_matrix->_maxModes];
   _modes = new qint8[_matrix->_sampleSize*_matrix->_maxModes];
}






CorrelationMatrix::Pair::Pair(const CorrelationMatrix* matrix):
   _cMatrix(matrix)
{
   _correlations = new Correlation[_matrix->_correlationSize*_matrix->_maxModes];
   _modes = new qint8[_matrix->_sampleSize*_matrix->_maxModes];
}






void CorrelationMatrix::Pair::write(int x, int y)
{
   if ( x < 0 || y < 0 || x >= y )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to write invalid pair %1,%2 with gene size of %3.").arg(x).arg(y)
                   .arg(_matrix->_geneSize));
      throw e;
   }
   _matrix->writeCorrelation(x,y,_correlations,_modes,_modeSize);
}






void CorrelationMatrix::Pair::read(int x, int y) const
{
   if ( x < 0 || y < 0 || x >= y )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to read invalid pair %1,%2 with gene size of %3.").arg(x).arg(y)
                   .arg(_cMatrix->_geneSize));
      throw e;
   }
   _cMatrix->readCorrelation(x,y,_correlations,_modes,_modeSize);
}






qint8 CorrelationMatrix::Pair::getModeSize() const
{
   return _modeSize;
}






CorrelationMatrix::Correlation& CorrelationMatrix::Pair::at(int mode, int index)
{
   if ( mode < 0 || mode > _modeSize || index < 0 || index > _matrix->_correlationSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Invalid mode of %1 and index of %2 with max modes of %3 and correlation "
                      "size of %4").arg(mode).arg(index).arg(_matrix->_maxModes)
                   .arg(_matrix->_correlationSize));
      throw e;
   }
   return _correlations[(mode*_matrix->_correlationSize)+index];
}






const CorrelationMatrix::Correlation& CorrelationMatrix::Pair::at(int mode, int index) const
{
   if ( mode < 0 || mode > _modeSize || index < 0 || index > _cMatrix->_correlationSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Invalid mode of %1 and index of %2 with max modes of %3 and correlation "
                      "size of %4").arg(mode).arg(index).arg(_cMatrix->_maxModes)
                   .arg(_cMatrix->_correlationSize));
      throw e;
   }
   return _correlations[(mode*_cMatrix->_correlationSize)+index];
}






qint8& CorrelationMatrix::Pair::mode(int mode, int index)
{
   if ( mode < 0 || mode > _modeSize || index < 0 || index > _matrix->_sampleSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Invalid mode of %1 and index of %2 with max modes of %3 and sample size of"
                      " %4").arg(mode).arg(index).arg(_matrix->_maxModes)
                   .arg(_matrix->_sampleSize));
      throw e;
   }
   return _modes[(mode*_matrix->_sampleSize)+index];
}






const qint8& CorrelationMatrix::Pair::mode(int mode, int index) const
{
   if ( mode < 0 || mode > _modeSize || index < 0 || index > _cMatrix->_sampleSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Invalid mode of %1 and index of %2 with max modes of %3 and sample size of"
                      " %4").arg(mode).arg(index).arg(_cMatrix->_maxModes)
                   .arg(_cMatrix->_sampleSize));
      throw e;
   }
   return _modes[(mode*_cMatrix->_sampleSize)+index];
}
