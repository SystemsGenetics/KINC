#include <ace/core/metadata.h>

#include "correlationmatrix.h"



using namespace std;






void CorrelationMatrix::readData()
{
   // seek to beginning of data and make sure it worked
   if ( !seek(0) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed calling seek() on data object file."));
      throw e;
   }

   // read in all header values
   stream() >> _geneSize >> _sampleSize >> _correlationSize >> _maxModes;

   // make sure reading worked
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
   // initialize gene, correlation, and mode size
   qint64 geneSize {_geneSize};
   qint64 correlationSize = _correlationSize*_maxModes*sizeof(Correlation);
   qint64 modeSize = _sampleSize*_maxModes*sizeof(qint8);

   // compute total size of data and return
   return (geneSize*(geneSize-1)/2)*(correlationSize+modeSize+sizeof(qint8));
}






void CorrelationMatrix::newData()
{
   // seek to beginning of data and make sure it worked
   if ( !seek(0) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed calling seek() on data object file."));
      throw e;
   }

   // write out all header information
   stream() << _geneSize << _sampleSize << _correlationSize << _maxModes;

   // make sure writing worked
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
   // check to see if pre-allocation is requested
   if ( preAllocate )
   {
      // seek to beginning of data and make sure it worked
      if ( !seek(0) )
      {
         E_MAKE_EXCEPTION(e);
         e.setTitle(tr("File IO Error"));
         e.setDetails(tr("Failed calling seek() on data object file."));
         throw e;
      }

      // allocate total size needed in data and make sure it worked
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
   // orientation is not used
   Q_UNUSED(orientation);

   // if role is not display return nothing
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // get map to metadata root and check if genes key exists
   const EMetadata::Map* map {meta().toObject()};
   if ( map->contains("genes") )
   {
      // get genes metadata and make sure it is an array
      EMetadata* genes {(*map)["genes"]};
      if ( genes->isArray() )
      {
         // make sure section is within limits of gene name array
         if ( section >= 0 && section < genes->toArray()->size() )
         {
            // return gene name
            return genes->toArray()->at(section)->toVariant();
         }
      }
   }

   // no gene found return nothing
   return QVariant();
}






int CorrelationMatrix::rowCount(const QModelIndex& parent) const
{
   // return gene size
   Q_UNUSED(parent);
   return _geneSize;
}






int CorrelationMatrix::columnCount(const QModelIndex& parent) const
{
   // return gene size
   Q_UNUSED(parent);
   return _geneSize;
}






QVariant CorrelationMatrix::data(const QModelIndex& index, int role) const
{
   // if role is not display return nothing
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // if row and column are equal return one
   if ( index.row() == index.column() )
   {
      return 1;
   }

   // get constant pair and read in values
   const Pair pair(this);
   int x {index.row()};
   int y {index.column()};
   if ( y > x )
   {
      swap(x,y);
   }
   pair.read(x,y);

   // construct pair information string
   QString ret;
   for (int i = 0; i < pair.getModeSize() ;++i)
   {
      // begin new mode line
      ret.append(tr("Mode %1: ").arg(i));
      for (int j = 0; j < _correlationSize ;)
      {
         // add new correlation
         ret.append(QString::number(pair.at(i,j)));

         // add comma if this is not last correlation
         if ( ++j < _correlationSize )
         {
            ret.append(", ");
         }
      }

      // finish mode line with end of line
      ret.append("\n");
   }

   // return pair information string
   return ret;
}






void CorrelationMatrix::initialize(const EMetadata& geneNames, qint32 sampleSize
                                   , qint8 correlationSize, qint8 maxModes)
{
   // make sure gene names metadata is an array
   if ( !geneNames.isArray() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Gene names metadata is not correct array type."));
      throw e;
   }

   // get map of metadata root and make copy of gene names
   EMetadata::Map* map {meta().toObject()};
   map->insert("genes",new EMetadata(geneNames));

   // save gene, sample, correlation sizes and max modes
   _geneSize = geneNames.toArray()->size();
   _sampleSize = sampleSize;
   _correlationSize = correlationSize;
   _maxModes = maxModes;
}






void CorrelationMatrix::increment(int &x, int &y)
{
   // increment y and check if it is equal or bigger than x
   if ( ++y >= x )
   {
      // increment x and reset y to zero
      ++x;
      y = 0;
   }
}






void CorrelationMatrix::readPair(int x, int y, Correlation* correlations, qint8* modes
                                 , qint8& modeSize) const
{
   // initialize index, correlation size, and mode size
   qint64 bX {x};
   qint64 bY {y};
   qint64 index {(bX*(bX-1)/2)+bY};
   qint64 correlationSize = _correlationSize*_maxModes*sizeof(Correlation);
   qint64 iModeSize = _sampleSize*_maxModes*sizeof(qint8);

   // seek to beginning of pair information and make sure it worked
   if ( !seek(DATA_OFFSET+(index*(correlationSize+iModeSize+sizeof(qint8)))) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed calling seek() on data object file."));
      throw e;
   }

   // read in pair information
   stream() >> modeSize;
   stream().read(correlations,_correlationSize);
   stream().read(modes,iModeSize);

   // make sure reading worked
   if ( !stream() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed reading from data object file."));
      throw e;
   }
}






void CorrelationMatrix::writePair(int x, int y, const Correlation* correlations, const qint8* modes
                                  , qint8 modeSize)
{
   // initialize index, correlation size, and mode size
   qint64 bX {x};
   qint64 bY {y};
   qint64 index {(bX*(bX-1)/2)+bY};
   qint64 correlationSize = _correlationSize*_maxModes*sizeof(Correlation);
   qint64 iModeSize = _sampleSize*_maxModes*sizeof(qint8);

   // seek to beginning of pair information and make sure it worked
   if ( !seek(DATA_OFFSET+(index*(correlationSize+iModeSize+sizeof(qint8)))) )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("File IO Error"));
      e.setDetails(tr("Failed calling seek() on data object file."));
      throw e;
   }

   // write out pair information
   stream() << modeSize;
   stream().write(correlations,_correlationSize);
   stream().write(modes,iModeSize);

   // make sure writing worked
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
   // make sure matrix pointer is valid
   if ( !matrix )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to create gene pair with null correlation matrix pointer."));
      throw e;
   }

   // create new arrays for correlation values and mode masks
   _correlations = new Correlation[_matrix->_correlationSize*_matrix->_maxModes];
   _modes = new qint8[_matrix->_sampleSize*_matrix->_maxModes];
}






CorrelationMatrix::Pair::Pair(const CorrelationMatrix* matrix):
   _cMatrix(matrix)
{
   // make sure matrix pointer is valid
   if ( !matrix )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to create gene pair with null correlation matrix pointer."));
      throw e;
   }

   // create new arrays for correlation values and mode masks
   _correlations = new Correlation[_cMatrix->_correlationSize*_cMatrix->_maxModes];
   _modes = new qint8[_cMatrix->_sampleSize*_cMatrix->_maxModes];
}






void CorrelationMatrix::Pair::write(int x, int y)
{
   // make sure xy gene pair is valid
   if ( x < 0 || y < 0 || y >= x )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to write invalid pair %1,%2 with gene size of %3.").arg(x).arg(y)
                   .arg(_matrix->_geneSize));
      throw e;
   }

   // write pair to matrix
   _matrix->writePair(x,y,_correlations,_modes,_modeSize);
}






void CorrelationMatrix::Pair::read(int x, int y) const
{
   // make sure xy gene pair is valid
   if ( x < 0 || y < 0 || y >= x )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Attempting to read invalid pair %1,%2 with gene size of %3.").arg(x).arg(y)
                   .arg(_cMatrix->_geneSize));
      throw e;
   }

   // read pair from matrix
   _cMatrix->readPair(x,y,_correlations,_modes,_modeSize);
}






CorrelationMatrix::Correlation& CorrelationMatrix::Pair::at(int mode, int index)
{
   // make sure mode and index are valid
   if ( mode < 0 || mode > _modeSize || index < 0 || index > _matrix->_correlationSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Invalid mode of %1 and index of %2 with max modes of %3 and correlation "
                      "size of %4").arg(mode).arg(index).arg(_matrix->_maxModes)
                   .arg(_matrix->_correlationSize));
      throw e;
   }

   // return correlation value
   return _correlations[(mode*_matrix->_correlationSize)+index];
}






const CorrelationMatrix::Correlation& CorrelationMatrix::Pair::at(int mode, int index) const
{
   // make sure mode and index are valid
   if ( mode < 0 || mode > _modeSize || index < 0 || index > _cMatrix->_correlationSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Invalid mode of %1 and index of %2 with max modes of %3 and correlation "
                      "size of %4").arg(mode).arg(index).arg(_cMatrix->_maxModes)
                   .arg(_cMatrix->_correlationSize));
      throw e;
   }

   // return correlation value
   return _correlations[(mode*_cMatrix->_correlationSize)+index];
}






qint8& CorrelationMatrix::Pair::mode(int mode, int index)
{
   // make sure mode and index are valid
   if ( mode < 0 || mode > _modeSize || index < 0 || index > _matrix->_sampleSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Invalid mode of %1 and index of %2 with max modes of %3 and sample size of"
                      " %4").arg(mode).arg(index).arg(_matrix->_maxModes)
                   .arg(_matrix->_sampleSize));
      throw e;
   }

   // return mode sample mask
   return _modes[(mode*_matrix->_sampleSize)+index];
}






const qint8& CorrelationMatrix::Pair::mode(int mode, int index) const
{
   // make sure mode and index are valid
   if ( mode < 0 || mode > _modeSize || index < 0 || index > _cMatrix->_sampleSize )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Invalid mode of %1 and index of %2 with max modes of %3 and sample size of"
                      " %4").arg(mode).arg(index).arg(_cMatrix->_maxModes)
                   .arg(_cMatrix->_sampleSize));
      throw e;
   }

   // return mode sample mask
   return _modes[(mode*_cMatrix->_sampleSize)+index];
}
