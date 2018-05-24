#include <ace/core/emetaarray.h>
#include <ace/core/emetaobject.h>

#include "ccmatrix.h"



using namespace std;
using namespace Pairwise;






void CCMatrix::Pair::addCluster(int amount) const
{
   // keep adding a new list of sample masks for given amount
   while ( amount-- > 0 )
   {
      _sampleMasks.append(QVector<qint8>(_cMatrix->_sampleSize, 0));
   }
}






QString CCMatrix::Pair::toString() const
{
   // if there are no clusters return empty string
   if ( _sampleMasks.isEmpty() )
   {
      return "";
   }

   // initialize list of strings and iterate through all clusters
   QStringList ret;
   for (const auto& sampleMask : _sampleMasks)
   {
      // initialize list of strings for sample mask and iterate through each sample
      QString clusterString("(");
      for (const auto& sample : sampleMask)
      {
         // add new sample token as hexadecimal allowing 16 different possible values
         switch (sample)
         {
         case 0:
         case 1:
         case 2:
         case 3:
         case 4:
         case 5:
         case 6:
         case 7:
         case 8:
         case 9:
            clusterString.append(QString::number(sample));
            break;
         case 10:
            clusterString.append("A");
            break;
         case 11:
            clusterString.append("B");
            break;
         case 12:
            clusterString.append("C");
            break;
         case 13:
            clusterString.append("D");
            break;
         case 14:
            clusterString.append("E");
            break;
         case 15:
            clusterString.append("F");
            break;
         }
      }

      // join all cluster string into one string
      ret << clusterString.append(')');
   }

   // join all clusters and return as string
   return ret.join(',');
}






void CCMatrix::Pair::writeCluster(EDataStream &stream, int cluster)
{
   // make sure cluster value is within range
   if ( cluster >= 0 && cluster < _sampleMasks.size() )
   {
      // write each sample to output stream
      auto& samples {_sampleMasks.at(cluster)};

      for ( int i = 0; i < samples.size(); i += 2 )
      {
         qint8 value {(qint8)(samples[i] & 0x0F)};

         if ( i + 1 < samples.size() )
         {
            value |= (samples[i + 1] << 4);
         }

         stream << value;
      }
   }
}






void CCMatrix::Pair::readCluster(const EDataStream &stream, int cluster) const
{
   // make sure cluster value is within range
   if ( cluster >= 0 && cluster < _sampleMasks.size() )
   {
      // read each sample from input stream
      auto& samples {_sampleMasks[cluster]};

      for ( int i = 0; i < samples.size(); i += 2 )
      {
         qint8 value;
         stream >> value;

         samples[i] = value & 0x0F;

         if ( i + 1 < samples.size() )
         {
            samples[i + 1] = (value >> 4) & 0x0F;
         }
      }
   }
}






QAbstractTableModel* CCMatrix::model()
{
   return nullptr;
}






QVariant CCMatrix::headerData(int section, Qt::Orientation orientation, int role) const
{
   // orientation is not used
   Q_UNUSED(orientation);

   // if role is not display return nothing
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // get genes metadata and make sure it is an array
   const EMetadata& genes {geneNames()};
   if ( genes.isArray() )
   {
      // make sure section is within limits of gene name array
      if ( section >= 0 && section < genes.toArray().size() )
      {
         // return gene name
         return genes.toArray().at(section).toString();
      }
   }

   // no gene found return nothing
   return QVariant();
}






int CCMatrix::rowCount(const QModelIndex&) const
{
   return geneSize();
}






int CCMatrix::columnCount(const QModelIndex&) const
{
   return geneSize();
}






QVariant CCMatrix::data(const QModelIndex &index, int role) const
{
   // if role is not display return nothing
   if ( role != Qt::DisplayRole )
   {
      return QVariant();
   }

   // if row and column are equal return empty string
   if ( index.row() == index.column() )
   {
      return "";
   }

   // get constant pair and read in values
   const Pair pair(this);
   int x {index.row()};
   int y {index.column()};
   if ( y > x )
   {
      swap(x,y);
   }
   pair.read({x,y});

   // Return value of gene pair as a string
   return pair.toString();
}






void CCMatrix::initialize(const EMetadata &geneNames, int maxClusterSize, const EMetadata &sampleNames)
{
   // make sure sample names is an array and is not empty
   if ( !sampleNames.isArray() || sampleNames.toArray().isEmpty() )
   {
      E_MAKE_EXCEPTION(e);
      e.setTitle(tr("Domain Error"));
      e.setDetails(tr("Sample names metadata is not an array or is empty."));
      throw e;
   }

   // save sample names to metadata
   // meta().toObject().insert("samples", sampleNames);

   // save sample size and initialize base class
   _sampleSize = sampleNames.toArray().size();
   Matrix::initialize(geneNames, maxClusterSize, (_sampleSize + 1) / 2 * sizeof(qint8), DATA_OFFSET);
}






EMetadata CCMatrix::sampleNames() const
{
   return meta().toObject().at("samples");
}
