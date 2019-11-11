#include "conditionspecificclustersmatrix_pair.h"
#include "pairwise_index.h"
#include <sstream>
#include <string>



/*!
 * Initialize a new CSM pair.
 *
 * @param matrix The parent matrix of the model.
 */
CSMatrix::Pair::Pair(CSMatrix* matrix) : Matrix::Pair(matrix), _cMatrix(matrix)
{
    EDEBUG_FUNC(this, matrix, testType);
}



/*!
 * Initialize a new CSM pair.
 *
 * @param matrix The parent matrix of the model.
 */
CSMatrix::Pair::Pair(const CSMatrix* matrix) : Matrix::Pair(matrix), _cMatrix(matrix)
{
    EDEBUG_FUNC(this, matrix, testType);
}



/*!
 * Clear all of the clusters out of a pair.
 */
void CSMatrix::Pair::clearClusters() const
{
    EDEBUG_FUNC(this);

    _pValues.clear();
    _r2.clear();
}



/*!
 * Add up to the max amount of clusters to the pair.
 *
 * @param amount The number of cluster you want to add.
 */
void CSMatrix::Pair::addCluster(int amount) const
{
    EDEBUG_FUNC(this, amount);

    while ( amount-- > 0 )
    {
        _pValues.append(QVector<double>());
        _r2.append(QVector<double>());
    }
}



/*!
 * Add up to the max amount of clusters to the pair. It also
 * resizes the vector to the chosen size, usually the test size in the case of CSM.
 *
 * @param amount The number of cluster you want to add.
 */
void CSMatrix::Pair::addCluster(int amount, int size) const
{
    EDEBUG_FUNC(this, amount, size);

    while ( amount-- > 0 )
    {
        _pValues.append(QVector<double>());
        _pValues.last().resize(size);
        _r2.append(QVector<double>());
        _r2.last().resize(size);
    }
}



/*!
 * Query for the size of a cluster.
 *
 * @return An integer representation of the size of a cluster.
 */
int CSMatrix::Pair::clusterSize() const
{
    EDEBUG_FUNC(this);

    return _pValues.size();
}



/*!
 * Query about a pairs contents.
 *
 * @return True if the pair is empty and false otherwise.
 */
bool CSMatrix::Pair::isEmpty() const
{
    EDEBUG_FUNC(this);

    return _pValues.isEmpty();
}



/*!
 * Convert the contents of the pair into a string.
 *
 * @return The string representation of the pair.
 */
QString CSMatrix::Pair::toString() const
{
    EDEBUG_FUNC(this);

    QString outputString;

    // if there is at least one cluster
    if ( !_pValues.isEmpty() )
    {
        // for each cluster
        for ( int i = 0; i < _pValues.size(); i++ )
        {
            outputString += "Cluster: ";
            outputString += QString::number(i + 1);
            outputString += "\n";

            for ( int j = 0; j < _pValues.at(i).size(); j++ )
            {
                double pvalue = _pValues.at(i).at(j);
                double r2 = _r2.at(i).at(j);

                outputString += "\t";
                outputString += _cMatrix->getTestName(j);
                outputString += ": p-value=";
                outputString += QString::number(pvalue);

                // If r2 is set to NaN then this was not a
                // test that generated the r2 calculation, so leave it out.
                if ( !qIsNaN(_r2.at(i).at(j)) )
                {
                    outputString += ", r2=";
                    outputString += QString::number(r2);
                }
                outputString += "\n";
            }
        }
    }

    return outputString;
}



/*!
 * Query for a value in a cluster in a sample.
 *
 * @param cluster The numeric index of the cluster
 *
 * @param gene The numeric index of the gene
 *
 * @param type The type of value to retrieve: 'pvalue' or 'r2'.
 *
 * @return The data at the quieried location.
 */
const double& CSMatrix::Pair::at(int cluster, int gene, QString type) const
{
    EDEBUG_FUNC(this, cluster, gene);

    if ( QString::compare(type, "pvalue") == 0 )
    {
        return _pValues.at(cluster).at(gene);
    }

    if ( QString::compare(type, "r2") == 0 )
    {
        return _r2.at(cluster).at(gene);
    }
    else
    {
        E_MAKE_EXCEPTION(e);
        e.setTitle(tr("Invalid type"));
        e.setDetails(tr("CSMatrix::Pair::at. The type argument must be a 0 or 1"));
        throw e;
    }
}



/*!
 * Query for a value in a cluster in a sample.
 *
 * @param cluster The numeric index of the cluster
 *
 * @param gene The numeric index of the gene
 *
 * @param type The type of value to retrieve: 'pvalue' or 'r2'.
 *
 * @return The data at the quieried location.
 */
double& CSMatrix::Pair::at(int cluster, int gene, QString type)
{
    EDEBUG_FUNC(this, cluster, gene);

    if ( QString::compare(type, "pvalue") == 0 )
    {
        return _pValues[cluster][gene];
    }

    if ( QString::compare(type, "r2") == 0 )
    {
        return _r2[cluster][gene];
    }
    else
    {
        E_MAKE_EXCEPTION(e);
        e.setTitle(tr("Invalid type"));
        e.setDetails(tr("CSMatrix::Pair::at. The type argument must be a 0 or 1"));
        throw e;
    }
}



/*!
 * Write cluster data into the file ACE sets up.
 *
 * @param stream The file stream to write into.
 *
 * @param cluster The cluster to write to the stream.
 */
void CSMatrix::Pair::writeCluster(EDataStream& stream, int cluster)
{
    EDEBUG_FUNC(this,&stream,cluster);

    // make sure cluster value is within range
    if ( cluster >= 0 && cluster < _pValues.size() )
    {
        // write each pvalue and r2 value to the output stream
        for ( qint32 i = 0; i < _cMatrix->_testcount; i ++ )
        {
            stream << _pValues[cluster][i] << _r2[cluster][i];
        }
    }
}



/*!
 * Read cluster data into the file ACE sets up.
 *
 * @param stream The file stream to read into.
 *
 * @param cluster The cluster to read from the stream.
 */
void CSMatrix::Pair::readCluster(const EDataStream& stream, int cluster) const
{
    EDEBUG_FUNC(this,&stream,cluster);

    // make sure cluster value is within range
    if ( cluster >= 0 && cluster < _pValues.size() )
    {
        // read each pvalue from input stream
        for ( int i = 0; i < _cMatrix->_testcount; i++ )
        {
            double pvalue = 0.0;
            double r2 = 0.0;
            stream >> pvalue >> r2;
            _pValues[cluster].append(pvalue);
            _r2[cluster].append(r2);
        }
    }
}
