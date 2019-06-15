#include "ccmatrix_pair.h"


/*!
 * Add one or more clusters to this pair.
 *
 * @param amount
 */
void CCMatrix::Pair::addCluster(int amount) const
{
   EDEBUG_FUNC(this,amount);

   // keep adding a new list of sample masks for given amount
   while ( amount-- > 0 )
   {
      _sampleMasks.append(QVector<qint8>(_cMatrix->_sampleSize, 0));
   }
}






/*!
 * Return the string representation of this pair, which is a comma-delimited
 * string of each sample mask in the pair.
 */
QString CCMatrix::Pair::toString() const
{
   EDEBUG_FUNC(this);

   // if there are no clusters return empty string
   if ( _sampleMasks.isEmpty() )
   {
      return "";
   }

   // initialize list of strings and iterate through all clusters
   QStringList ret;
   for ( const auto& sampleMask : _sampleMasks )
   {
      // initialize list of strings for sample mask and iterate through each sample
      QString clusterString("(");
      for ( const auto& sample : sampleMask )
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






/*!
 * Write a cluster in the iterator's pairwise data to the data object file.
 *
 * @param stream
 * @param cluster
 */
void CCMatrix::Pair::writeCluster(EDataStream& stream, int cluster)
{
   EDEBUG_FUNC(this,&stream,cluster);

   // make sure cluster value is within range
   if ( cluster >= 0 && cluster < _sampleMasks.size() )
   {
      // write each sample to output stream
      auto& samples {_sampleMasks.at(cluster)};

      for ( int i = 0; i < samples.size(); i += 2 )
      {
         qint8 value {static_cast<qint8>(samples[i] & 0x0F)};

         if ( i + 1 < samples.size() )
         {
            value |= static_cast<qint8>(samples[i + 1] << 4);
         }

         stream << value;
      }
   }
}






/*!
 * Read a cluster from the data object file into memory.
 *
 * @param stream
 * @param cluster
 */
void CCMatrix::Pair::readCluster(const EDataStream& stream, int cluster) const
{
   EDEBUG_FUNC(this,&stream,cluster);

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
