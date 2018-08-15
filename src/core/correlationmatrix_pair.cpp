#include "correlationmatrix_pair.h"



void CorrelationMatrix::Pair::addCluster(int amount) const
{
   // keep adding a new list of floats for given amount
   while ( amount-- > 0 )
   {
      _correlations.append(QVector<float>(_cMatrix->_correlationSize, NAN));
   }
}






QString CorrelationMatrix::Pair::toString() const
{
   // if there are no correlations simply return null
   if ( _correlations.isEmpty() )
   {
      return tr("");
   }

   // initialize list of strings and iterate through all clusters
   QStringList ret;
   for (const auto& cluster : _correlations)
   {
      // initialize list of strings for cluster and iterate through each correlation
      QStringList clusterStrings;
      for (const auto& correlation : cluster)
      {
         // add correlation value as string
         clusterStrings << QString::number(correlation);
      }

      // join all cluster strings into one string
      ret << clusterStrings.join(',');
   }

   // join all clusters and return as string
   return ret.join(',');
}






void CorrelationMatrix::Pair::writeCluster(EDataStream& stream, int cluster)
{
   // make sure cluster value is within range
   if ( cluster >= 0 && cluster < _correlations.size() )
   {
      // write correlations per cluster to output stream
      for (const auto& correlation : _correlations.at(cluster))
      {
         stream << correlation;
      }
   }
}






void CorrelationMatrix::Pair::readCluster(const EDataStream& stream, int cluster) const
{
   // make sure cluster value is within range
   if ( cluster >= 0 && cluster < _correlations.size() )
   {
      // read correlations per cluster from input stream
      for (int i = 0; i < _cMatrix->_correlationSize ;++i)
      {
         float value;
         stream >> value;
         _correlations[cluster][i] = value;
      }
   }
}
