#include <math.h>
#include "extract_networkwriter.h"
#include "expressionmatrix_gene.h"



/*!
 * The constructor for the NetworkWriter class.
 *
 * @param emx  The expression matrix.
 * @param cmx  The correlation matrix.
 * @param ccm  The cluster composition matrix.
 * @param csm  The condition-specific matrix.
 * @param output  The output network file.
 */
Extract::NetworkWriter::NetworkWriter(
    ExpressionMatrix* emx,
    CorrelationMatrix* cmx,
    CCMatrix* ccm,
    CSMatrix* csm,
    QFile* output):
    _emx(emx),
    _cmx(cmx),
    _cmxPair(cmx),
    _ccm(ccm),
    _ccmPair(ccm),
    _csm(csm),
    _csmPair(csm),
    _output(output)
{
    EDEBUG_FUNC(this,stream,emx,cmx,ccm,csm);

    // initialize output file stream
    _stream.setDevice(_output);
    _stream.setRealNumberPrecision(8);
}



/*!
 * Read the next gene pair and return the cluster size.
 */
int Extract::NetworkWriter::readNext()
{
    EDEBUG_FUNC(this,cmx_index);

    // read the next pair in the cmx
    _cmxPair.readNext();

    // read the next pair in the ccm (should always match cmx)
    if ( _ccm )
    {
        _ccmPair.readNext();

        if ( _cmxPair.index() != _ccmPair.index() )
        {
            qInfo() << "warning: cmx and ccm are out of sync at cmx coordinate ("
                    << _cmxPair.index().getX() << "," << _cmxPair.index().getY() <<").";
        }
    }

    // read the next pair in the csm (should always match cmx)
    if ( _csm )
    {
        _csmPair.readNext();

        if ( _cmxPair.index() != _csmPair.index() )
        {
            qInfo() << "warning: cmx and ccm are out of sync at cmx coordinate ("
                    << _cmxPair.index().getX() << "," << _cmxPair.index().getY() <<").";
        }
    }

    return _cmxPair.clusterSize();
}



/*!
 * Returns the name of the first gene in a gene pair. The gene pair must first
 * be set using the setPair() function.
 */
QString Extract::NetworkWriter::getEdgeGene1() const
{
    EDEBUG_FUNC(this);

    return _cmx->geneNames().at(_cmxPair.index().getX()).toString();
}



/*!
 * Returns the name of the second gene in a gene pair. The gene pair must first
 * be set using the setPair() function.
 */
QString Extract::NetworkWriter::getEdgeGene2() const
{
    EDEBUG_FUNC(this);

    return _cmx->geneNames().at(_cmxPair.index().getY()).toString();
}



/*!
 * Returns the similarity score of the cluster specified by
 * the cluster index.  The cluster must belong to the current gene pair.
 * The gene pair must first be set using the setPair() function.
 *
 * @param clusterIndex  The index of the cluster in the gene pair.
 */
float Extract::NetworkWriter::getEdgeSimilarity(int clusterIndex) const
{
    EDEBUG_FUNC(this,clusterIndex);

    return _cmxPair.at(clusterIndex);
}



/*!
 * Retrieves the sample string of the specified cluster
 * in the current gene pair. The gene pair must first
 * be set using the setPair() function.
 *
 * @param clusterIndex The index of the cluster in the current gene pair.
 *
 * \return  The sample string of characters indicating the role
 * of the sample in the cluster.
 */
QString Extract::NetworkWriter::getEdgeSampleString(int clusterIndex) const
{
    EDEBUG_FUNC(this,clusterIndex);

    // Initialize the sample string with zeros.
    QString sampleMask(_ccm->sampleSize(), '0');

    // If cluster data exists then use it.
    if ( _ccmPair.clusterSize() > 0 )
    {
        // Write sample mask to string.
        for ( int i = 0; i < _ccm->sampleSize(); i++ )
        {
            sampleMask[i] = '0' + _ccmPair.at(clusterIndex, i);
        }
    }

    // Otherwise use expression data if provided.
    else if ( _emx )
    {
        // Read in gene expressions.
        ExpressionMatrix::Gene gene1(_emx);
        ExpressionMatrix::Gene gene2(_emx);

        gene1.read(_cmxPair.index().getX());
        gene2.read(_cmxPair.index().getY());

        // Determine sample mask, summary statistics from expression data.
        for ( int i = 0; i < _emx->sampleSize(); ++i )
        {
            if ( isnan(gene1.at(i)) || isnan(gene2.at(i)) )
            {
                sampleMask[i] = '9';
            }
            else
            {
                sampleMask[i] = '1';
            }
        }
    }

    // Otherwise throw an error.
    else
    {
        E_MAKE_EXCEPTION(e)
        e.setTitle(QObject::tr("Invalid Input"));
        e.setDetails(QObject::tr("Expression Matrix was not provided but Cluster Matrix is missing sample data."));
        throw e;
    }

    return sampleMask;
}



/*!
 * Retrieves the total number of samples in the specified cluster
 * in the current gene pair.
 *
 * @param sampleMask  The sample string as returned by getEdgeSampleString().
 */
int Extract::NetworkWriter::getEdgeNumSamples(QString sampleMask) const
{
    EDEBUG_FUNC(this,clusterIndex);

    return sampleMask.count("1");
}



/*!
 * If the condition-specific matrix (csm) is present then this
 * function gets the list of tests that were performed and stores
 * them in a variable local to this class for easy iteration later.
 * This function should be called in the initialize() function of
 * a child NetworkWriter implementation if the test data will be written
 * in the output file.
 */
void Extract::NetworkWriter::setTestNames()
{
    EDEBUG_FUNC(this);

    if (!_csm)
    {
        return;
    }

    for ( int i = 0; i < _csm->getTestCount(); i++ )
    {
        QString testName = _csm->getTestName(i);
        QString testType = _csm->getTestType(i);

        // There is always a p-value.
        _testNames.append(testName + + "_pVal");

        // If linear regression was performed then there should
        // be an R-squared value as well.
        if ( QString::compare(testType, "Quantitative", Qt::CaseInsensitive) == 0 ||
             QString::compare(testType, "Ordinal", Qt::CaseInsensitive) == 0 )
        {
            _testNames.append(testName + + "_RSqr");
        }
    }
}



/*!
 * Retrieves the test value for a specified test, from the specified
 * cluster of the current gene pair. The gene pair must first
 * be set using the setPair() function.
 *
 * @param clusterIndex  The index of the cluster in the current gene pair.
 * @param fullTestName The name of the test to retrieve the value for.
 */
float Extract::NetworkWriter::getEdgeTestValue(int clusterIndex, const QString& fullTestName) const
{
    EDEBUG_FUNC(this,clusterIndex,fullTestName);

    for ( int i = 0; i < _csm->getTestCount(); i++ )
    {
        QString testName = _csm->getTestName(i);

        if (fullTestName.contains("_pVal") && fullTestName.compare(testName + "_pVal") == 0)
        {
            return static_cast<float>(_csmPair.at(clusterIndex, i, "pvalue"));
        }

        if (fullTestName.contains("_RSqr") && fullTestName.compare(testName + "_RSqr") == 0)
        {
            return static_cast<float>(_csmPair.at(clusterIndex, i, "r2"));
        }
    }

    return static_cast<float>(qQNaN());
}



/*!
 * Checks to make sure that the output stream is ok.
 */
void Extract::NetworkWriter::checkStatus() const
{
    if ( _stream.status() != QTextStream::Ok )
    {
        E_MAKE_EXCEPTION(e)
        e.setTitle(QObject::tr("File IO Error"));
        e.setDetails(QObject::tr("Qt Text Stream encountered an unknown error."));
        throw e;
    }
}



/*!
 * The initialization function for the full text output. This
 * function adds the tab-delimited header to the output file.
 */
void Extract::FullNetworkWriter::initialize()
{
    EDEBUG_FUNC(this);

    setTestNames();

    // Start by writing the tab-delimited header.
    _stream
        << "Source"
        << "\t" << "Target"
        << "\t" << "Similarity_Score"
        << "\t" << "Interaction"
        << "\t" << "Cluster_Index"
        << "\t" << "Cluster_Size"
        << "\t" << "Samples";

    // Add each conditional test to the network file too.
    for (auto& testName : _testNames)
    {
        _stream << "\t" << testName;
    }

    _stream
        << "\n";
}



/*!
 * Writes a single edge corresponding to the given cluster of the current
 * gene pair to the full tab-delimited text output file. The gene pair
 * must first be set using the setPair() function.
 *
 * @param clusterIndex The index of cluster in the current gene pair.
 * @param passed The set of tests that passed. (unused)
 */
void Extract::FullNetworkWriter::writeEdgeCluster(int clusterIndex, QVector<QString>)
{
    EDEBUG_FUNC(this,clusterIndex);

    QString source = getEdgeGene1();
    QString target = getEdgeGene2();
    float correlation = getEdgeSimilarity(clusterIndex);
    QString interaction {"co"};
    int cluster_num = clusterIndex + 1;
    QString sampleMask = getEdgeSampleString(clusterIndex);
    int numSamples = getEdgeNumSamples(sampleMask);

    _stream
        << source
        << "\t" << target
        << "\t" << correlation
        << "\t" << interaction
        << "\t" << cluster_num
        << "\t" << numSamples
        << "\t" << sampleMask;

    // Add each conditional test to the network file too.
    for (auto& testName : _testNames)
    {
        float testValue = getEdgeTestValue(clusterIndex, testName);
        _stream << "\t" << testValue;
    }

    _stream
        << "\n";
}



/*!
 * The initialization function for the full text output. This
 * function adds the tab-delimited header to the output file.
 */
void Extract::TidyNetworkWriter::initialize()
{
    EDEBUG_FUNC(this);

    setTestNames();

    // Start by writing the tab-delimited header.
    _stream
        << "Source"
        << "\t" << "Target"
        << "\t" << "Similarity_Score"
        << "\t" << "Interaction"
        << "\t" << "Cluster_Index"
        << "\t" << "Cluster_Size"
        << "\t" << "Samples";

    if (_testNames.size() > 0)
    {
        _stream
            << "\t" << "Test_Name"
            << "\t" << "p_value"
            << "\t" << "r_squared";
    }

    _stream
        << "\n";
}



/*!
 * Writes a single edge corresponding to the given cluster of the current
 * gene pair to the full tab-delimited text output file. The gene pair
 * must first be set using the setPair() function.
 *
 * @param clusterIndex The index of cluster in the current gene pair.
 * @param passed the set of tests that passed.
 */
void Extract::TidyNetworkWriter::writeEdgeCluster(int clusterIndex, QVector<QString> passed)
{
    EDEBUG_FUNC(this,clusterIndex,passed);

    QString source = getEdgeGene1();
    QString target = getEdgeGene2();
    float correlation = getEdgeSimilarity(clusterIndex);
    QString interaction {"co"};
    int cluster_num = clusterIndex + 1;
    QString sampleMask = getEdgeSampleString(clusterIndex);
    int numSamples = getEdgeNumSamples(sampleMask);

    // If there are no tests then we are here becaues the cluster
    // passed the correlation threshold test, so write it out.
    if (_testNames.size() == 0)
    {
        _stream
            << source
            << "\t" << target
            << "\t" << correlation
            << "\t" << interaction
            << "\t" << cluster_num
            << "\t" << numSamples
            << "\t" << sampleMask
            << "\n";
    }

    // Otherwise, write out an edge for each passed test.
    else
    {
        for (auto& testName : passed)
        {
            int pVali = _testNames.lastIndexOf(testName + "_pVal");
            int rSqri = _testNames.lastIndexOf(testName + "_RSqr");

            float pVal = (pVali != -1)
                ? getEdgeTestValue(clusterIndex, _testNames[pVali])
                : qQNaN();

            float rSqr = (rSqri != -1)
                ? getEdgeTestValue(clusterIndex, _testNames[rSqri])
                : qQNaN();

            _stream
                << source
                << "\t" << target
                << "\t" << correlation
                << "\t" << interaction
                << "\t" << cluster_num
                << "\t" << numSamples
                << "\t" << sampleMask
                << "\t" << testName
                << "\t" << pVal
                << "\t" << rSqr
                << "\n";
        }
    }
}



/*!
 * The initialization function for the minimial text output. This
 * function adds the tab-delimited header to the output file.
 */
void Extract::MinimalNetworkWriter::initialize()
{
    EDEBUG_FUNC(this);

    // Start by writing the tab-delimited header.
    _stream
        << "Source"
        << "\t" << "Target"
        << "\t" << "Similarity_Score"
        << "\t" << "Cluster_Index"
        << "\t" << "Num_Clusters"
        << "\n";
}



/*!
 * Writes a single edge corresponding to the given cluster of the current
 * gene pair to the minimal tab-delimited text output file. The gene pair
 * must first be set using the setPair() function.
 *
 * @param clusterIndex The index of cluster in the current gene pair.
 * @param passed The set of tests that passed. (unused)
 */
void Extract::MinimalNetworkWriter::writeEdgeCluster(int clusterIndex, QVector<QString>)
{
    EDEBUG_FUNC(this,clusterIndex);

    QString source = getEdgeGene1();
    QString target = getEdgeGene2();
    float correlation = getEdgeSimilarity(clusterIndex);
    int cluster_num = clusterIndex + 1;

    _stream
        << source
        << "\t" << target
        << "\t" << correlation
        << "\t" << cluster_num
        << "\t" << _cmxPair.clusterSize()
        << "\n";
}



/*!
 * The initialization function for the GraphML output. This
 * function adds XML headers, and attribute keys to the
 * output file.
 */
void Extract::GMLNetworkWriter::initialize()
{
    EDEBUG_FUNC(this);

    setTestNames();

    // Write the XML header to file.
    _stream
            << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
            << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n"
            << "         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
            << "         xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n";

    // Add the attributes definitions
    _stream
            << "  <key id=\"sc\" for=\"edge\" attr.name=\"Similarity_Score\" attr.type=\"string\"/>\n"
            << "  <key id=\"int\" for=\"edge\" attr.name=\"Interaction\" attr.type=\"string\"/>\n"
            << "  <key id=\"ci\" for=\"edge\" attr.name=\"Cluster_Index\" attr.type=\"int\"/>\n"
            << "  <key id=\"cs\" for=\"edge\" attr.name=\"Cluster_Size\" attr.type=\"int\"/>\n"
            << "  <key id=\"s\" for=\"edge\" attr.name=\"Samples\" attr.type=\"string\"/>\n";

    // Add attributes for the tests.
    for (int i = 0; i < _testNames.size(); ++i)
    {
        _stream
            << "  <key id=\"t" << i << "\" for=\"node\" attr.name=\"" << _testNames[i] << "\" attr.type=\"double\"/>\n";
    }

    // Start the graph.
    _stream
            << "  <graph id=\"G01\" edgedefault=\"undirected\">\n";
}



/*!
 * Writes a single edge corresponding to the given cluster of the current
 * gene pair to the GraphML output file. The gene pair
 * must first be set using the setPair() function.
 *
 * @param clusterIndex The index of cluster in the current gene pair.
 * @param passed The set of tests that passed. (unused)
 */
void Extract::GMLNetworkWriter::writeEdgeCluster(int clusterIndex, QVector<QString>)
{
    EDEBUG_FUNC(this,clusterIndex);

    QString source = getEdgeGene1();
    QString target = getEdgeGene2();
    float correlation = getEdgeSimilarity(clusterIndex);
    QString interaction {"co"};
    int cluster_num = clusterIndex + 1;
    QString sampleMask = getEdgeSampleString(clusterIndex);
    int numSamples = getEdgeNumSamples(sampleMask);

    // If we've never seen this node before then add it to
    // our list and to the output file.
    if (!_nodes.contains(source))
    {
        _nodes.insert(source, true);

        _stream
            << "    <node id=\"" << source << "\"/>\n";
    }
    if (!_nodes.contains(target))
    {
        _nodes.insert(target, true);

        _stream
            << "    <node id=\"" << target << "\"/>\n";
    }

    // Start the edge and add in the attributes.
    _stream
            << "    <edge source=\"" << source << "\" target=\"" << target << "\">\n"
            << "      <data key=\"sc\">" << correlation << "</data>\n"
            << "      <data key=\"int\">co</data>\n"
            << "      <data key=\"ci\">" << cluster_num << "</data>\n"
            << "      <data key=\"cs\">" << numSamples << "</data>\n"
            << "      <data key=\"s\">" << sampleMask << "</data>\n";

    // Add each conditional test to the network file too.
    for (int i = 0; i < _testNames.size(); ++i)
    {
        float testValue = getEdgeTestValue(clusterIndex, _testNames[i]);

        _stream
            << "      <data key=\"t" << i << "\">" << testValue << "</data>\n";
    }

    // Finish out the edge.
    _stream
            << "    </edge>\n";
}



/*!
 * Terminates the GraphML output file.
 */
void Extract::GMLNetworkWriter::finish()
{
    EDEBUG_FUNC(this);

    // Write the XML header to file.
    _stream
            << "  </graph>\n"
            << "</graphml>\n";
}
