#include "conditionspecificclustersmatrix_pair.h"
#include "pairwise_index.h"
#include <sstream>
#include <string>
//





/*!
*  Implements an interface to initialize a new CSM pair.
*
* @param matrix The parent matrix of the model.
*/
CSMatrix::Pair::Pair(CSMatrix* matrix) : Matrix::Pair(matrix), _cMatrix(matrix)
{
    EDEBUG_FUNC(this,matrix);
}





/*!
*  Implements an interface to initialize a new CSM pair.
*
* @param matrix The parent matrix of the model.
*/
CSMatrix::Pair::Pair(const CSMatrix* matrix) : Matrix::Pair(matrix), _cMatrix(matrix)
{
    EDEBUG_FUNC(this,matrix);
}






/*!
*  Implements an interface to clear all of the clusters out of a pair.
*/
void CSMatrix::Pair::clearClusters() const
{
    EDEBUG_FUNC(this);
    _pValues.clear();
}





/*!
*  Implements an interface to add up to the max amount of clusters to the pair.
*
* @param amount The number of cluster you want to add.
*/
void CSMatrix::Pair::addCluster(int amount) const
{
    EDEBUG_FUNC(this, amount);
    while ( amount-- > 0 )
    {
       _pValues.append(QVector<double>());
    }
}






/*!
*  Implements an interface to add up to the max amount of clusters to the pair. It also
*  resizes the vector to the chosen size, usually the test size in the case of CSM.
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
    }
}





/*!
*  Implements an interface to quiery for the size of a cluster.
*
* @return An integer representation of the size of a cluster.
*/
int CSMatrix::Pair::clusterSize() const
{
    EDEBUG_FUNC(this);
    return  _pValues.size();
}





/*!
*  Implements an interface to quiery about a pairs contents.
*
* @return True if the pair is empty and false otherwise.
*/
bool CSMatrix::Pair::isEmpty() const
{
    EDEBUG_FUNC(this);
    return _pValues.isEmpty();
}





/*!
*  Implements an interface to convert the contents of the pair into a string.
*
* @param anxData The information about the annotation matrix, to help label
*        the string.
*
* @return The string representation of the pair.
*/
QString CSMatrix::Pair::toString() const
{
    EDEBUG_FUNC(this);

    QString outputString;

    //if there is at least one cluster
    if ( !_pValues.isEmpty() )
    {
        //for each cluster
        for ( int i = 0; i < _pValues.size(); i++ )
        {
            outputString+= "Cluster: ";
            outputString+= QString::number(i + 1);
            outputString+= "\n";
            for ( int j = 0; j < _pValues.at(i).size(); j++ )
            {
                outputString+= "\t";
                outputString+= _cMatrix->getTestName(j);
                outputString+= ": ";
                outputString+= QString::number(_pValues.at(i).at(j));
                outputString+= "\n";
            }
        }
    }
    return outputString;
}





/*!
*  Implements an interface to quiery for a value in a cluster in a sample.
*
* @return The data at the quieried location.
*/
const double& CSMatrix::Pair::at(int cluster, int gene) const
{
    EDEBUG_FUNC(this, cluster, gene);
    return _pValues.at(cluster).at(gene);
}





/*!
*  Implements an interface to quiery for a value in a cluster in a sample.
*
* @return The data at the quieried location.
*/
double& CSMatrix::Pair::at(int cluster, int gene)
{
    EDEBUG_FUNC(this, cluster, gene);
    return _pValues[cluster][gene];
}





/*!
*  Implements an interface to write cluster data into the file ACE sets up.
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
       // write each pvalue to output stream
       for ( qint32 i = 0; i < _cMatrix->_testcount; i ++ )
       {
           stream << _pValues[cluster][i];
       }
    }
}





/*!
*  Implements an interface to read cluster data into the file ACE sets up.
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
          double num = 0.0;
          stream >> num;
          _pValues[cluster].append(num);
       }
    }
}
