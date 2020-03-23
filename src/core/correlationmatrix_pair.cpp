#include "correlationmatrix_pair.h"



/*!
 * Add one or more clusters to this pair.
 *
 * @param amount
 */
void CorrelationMatrix::Pair::addCluster(int amount) const
{
    EDEBUG_FUNC(this,amount);

    // keep adding a new float for given amount
    while ( amount-- > 0 )
    {
        _correlations.append(NAN);
    }
}



/*!
 * Return the string representation of this pair, which is a comma-delimited
 * string of each correlation in the pair.
 */
QString CorrelationMatrix::Pair::toString() const
{
    EDEBUG_FUNC(this);

    // if there are no correlations simply return null
    if ( _correlations.isEmpty() )
    {
        return tr("");
    }

    // initialize list of strings and iterate through all clusters
    QStringList ret;
    for ( const auto& correlation : _correlations )
    {
        // add correlation value as string
        ret << QString::number(correlation);
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
void CorrelationMatrix::Pair::writeCluster(EDataStream& stream, int cluster)
{
    EDEBUG_FUNC(this,&stream,cluster);

    // make sure cluster value is within range
    if ( cluster >= 0 && cluster < _correlations.size() )
    {
        // write correlation to output stream
        stream << _correlations.at(cluster);
    }
}



/*!
 * Read a cluster from the data object file into memory.
 *
 * @param stream
 * @param cluster
 */
void CorrelationMatrix::Pair::readCluster(const EDataStream& stream, int cluster) const
{
    EDEBUG_FUNC(this,&stream,cluster);

    // make sure cluster value is within range
    if ( cluster >= 0 && cluster < _correlations.size() )
    {
        // read correlation from input stream
        stream >> _correlations[cluster];
    }
}
