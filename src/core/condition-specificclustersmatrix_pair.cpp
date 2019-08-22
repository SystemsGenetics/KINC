#include "condition-specificclustersmatrix_pair.h"
#include "pairwise_index.h"
#include <sstream>
#include <string>
//





/*!
*  Implements an interface to initialize a new CSCM pair.
*
* @param matrix The parent matrix of the model.
*/
CSCM::Pair::Pair(CSCM* matrix) : Matrix::Pair(matrix), _cMatrix(matrix)
{
    EDEBUG_FUNC(this,matrix);
}





/*!
*  Implements an interface to initialize a new CSCM pair.
*
* @param matrix The parent matrix of the model.
*/
CSCM::Pair::Pair(const CSCM* matrix) : Matrix::Pair(matrix), _cMatrix(matrix)
{
    EDEBUG_FUNC(this,matrix);
}






/*!
*  Implements an interface to clear all of the clusters out of a pair.
*/
void CSCM::Pair::clearClusters() const
{
    EDEBUG_FUNC(this);
    _pValues.clear();
}





/*!
*  Implements an interface to add up to the max amount of clusters to the pair.
*
* @param amount The number of cluster you want to add.
*/
void CSCM::Pair::addCluster(int amount) const
{
    EDEBUG_FUNC(this, amount);
    while ( amount-- > 0 )
    {
       _pValues.append(QVector<double>());
    }
}






/*!
*  Implements an interface to add up to the max amount of clusters to the pair. It also
*  resizes the vector to the chosen size, usually the test size in the case of CSCM.
*
* @param amount The number of cluster you want to add.
*/
void CSCM::Pair::addCluster(int amount, int size) const
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
int CSCM::Pair::clusterSize() const
{
    EDEBUG_FUNC(this);
    return  _pValues.size();
}





/*!
*  Implements an interface to quiery about a pairs contents.
*
* @return True if the pair is empty and false otherwise.
*/
bool CSCM::Pair::isEmpty() const
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
QString CSCM::Pair::toString() const
{
    EDEBUG_FUNC(this);

    std::string outputString;

    //if there is at least one cluster
    if(!_pValues.isEmpty())
    {
        //for each cluster
        for(int i = 0; i < _pValues.size(); i++)
        {
            outputString+= "Cluster: ";
            outputString+= static_cast<std::ostringstream*>( &(std::ostringstream() << i + 1) )->str();
            outputString+= "\n";
            for(int j = 0; j < _pValues.at(i).size(); j++)
            {
                outputString+= "\t";
                outputString+= _cMatrix->getTestName(j).toStdString();
                outputString+= ": ";
                outputString+= static_cast<std::ostringstream*>( &(std::ostringstream() << _pValues.at(i).at(j)) )->str();
                outputString+= "\n";
            }
        }
    }

    QString returnString;
    for(uint i = 0; i < outputString.size(); i++)
    {
        returnString.append(outputString.at(i));
    }
    return returnString;
}





/*!
*  Implements an interface to quiery for a value in a cluster in a sample.
*
* @return The data at the quieried location.
*/
const double& CSCM::Pair::at(int cluster, int gene) const
{
    EDEBUG_FUNC(this, cluster, gene);
    return _pValues.at(cluster).at(gene);
}





/*!
*  Implements an interface to quiery for a value in a cluster in a sample.
*
* @return The data at the quieried location.
*/
double& CSCM::Pair::at(int cluster, int gene)
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
void CSCM::Pair::writeCluster(EDataStream& stream, int cluster)
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
void CSCM::Pair::readCluster(const EDataStream& stream, int cluster) const
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
