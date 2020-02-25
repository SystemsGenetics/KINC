#include "cpmatrix_pair.h"



/*!
 * Add one or more clusters to this pair.
 *
 * @param amount
 */
void CPMatrix::Pair::addCluster(int amount) const
{
    EDEBUG_FUNC(this,amount);

    // keep adding a new component for given amount
    while ( amount-- > 0 )
    {
        _components.append(Component());
    }
}



/*!
 * Return the string representation of this pair, which is a comma-delimited
 * string of each component in the pair.
 */
QString CPMatrix::Pair::toString() const
{
    EDEBUG_FUNC(this);

    // if there are no components simply return null
    if ( _components.isEmpty() )
    {
        return tr("");
    }

    // initialize list of strings and iterate through all clusters
    QStringList ret;
    for ( const auto& component : _components )
    {
        // add component value as string
        ret << QString::asprintf("pi = %.3f mu = [ %.3f, %.3f ] sigma = [ %.3f, %.3f ; %.3f, %.3f ]",
            component.pi,
            component.mu.s[0],
            component.mu.s[1],
            component.sigma.s[0],
            component.sigma.s[1],
            component.sigma.s[2],
            component.sigma.s[3]
        );
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
void CPMatrix::Pair::writeCluster(EDataStream& stream, int cluster)
{
    EDEBUG_FUNC(this,&stream,cluster);

    // make sure cluster value is within range
    if ( cluster >= 0 && cluster < _components.size() )
    {
        // write component to output stream
        auto& component = _components.at(cluster);

        stream << component.pi;
        stream << component.mu.s[0];
        stream << component.mu.s[1];
        stream << component.sigma.s[0];
        stream << component.sigma.s[1];
        stream << component.sigma.s[2];
        stream << component.sigma.s[3];
    }
}



/*!
 * Read a cluster from the data object file into memory.
 *
 * @param stream
 * @param cluster
 */
void CPMatrix::Pair::readCluster(const EDataStream& stream, int cluster) const
{
    EDEBUG_FUNC(this,&stream,cluster);

    // make sure cluster value is within range
    if ( cluster >= 0 && cluster < _components.size() )
    {
        // read component from input stream
        auto& component = _components[cluster];

        stream >> component.pi;
        stream >> component.mu.s[0];
        stream >> component.mu.s[1];
        stream >> component.sigma.s[0];
        stream >> component.sigma.s[1];
        stream >> component.sigma.s[2];
        stream >> component.sigma.s[3];
    }
}
