#include "networkwriter.h"
#include "expressionmatrix_gene.h"
#include <math.h>



/*!
 * The constructor for the NetworkWriter class.
 *
 * @param emx  The expression matrix.
 * @param cmx  The correlation matrix.
 * @param ccm  The cluster composition matrix.
 * @param csm  The condition-specific matrix.
 * @param output  The output network file.
 */
NetworkWriter::NetworkWriter(
    ExpressionMatrix* emx,
    CorrelationMatrix* cmx,
    CCMatrix* ccm,
    CSMatrix* csm,
    QFile* output):
    _emx(emx),
    _cmx(cmx),
    _cmxPair(CorrelationMatrix::Pair(_cmx)),
    _ccm(ccm),
    _ccmPair(CCMatrix::Pair(ccm)),
    _csm(csm),
    _csmPair(CSMatrix::Pair(csm)),
    _output(output)
{
    EDEBUG_FUNC(this,stream,emx,cmx,ccm,csm);

    // initialize output file stream
    _stream.setDevice(_output);
    _stream.setRealNumberPrecision(8);
}



/*!
 * Returns the name of the first gene in a gene pair. The gene pair must first
 * be set using the setPair() function.
 */
QString NetworkWriter::getEdgeGene1() const
{
    EDEBUG_FUNC(this);

    EMetaArray geneNames {_cmx->geneNames()};
    return geneNames.at(_cmxPair.index().getX()).toString();
}



/*!
 * Returns the name of the second gene in a gene pair. The gene pair must first
 * be set using the setPair() function.
 */
QString NetworkWriter::getEdgeGene2() const
{
    EDEBUG_FUNC(this);

    EMetaArray geneNames {_cmx->geneNames()};
    return geneNames.at(_cmxPair.index().getY()).toString();
}



/*!
 * Returns the similarity score of the cluster specified by
 * the cluster index.  The cluster must belong to the current gene pair.
 * The gene pair must first be set using the setPair() function.
 *
 * @param cluster_index  The index of the cluster in the gene pair.
 */
float NetworkWriter::getEdgeSimilarity(int cluster_index) const
{
    EDEBUG_FUNC(this,cluster_index);

    return _cmxPair.at(cluster_index);
}



/*!
 * Read the next gene pair and return the cluster size.
 */
int NetworkWriter::readNext()
{
    EDEBUG_FUNC(this,cmx_index);

    _cmxPair.readNext();
    if ( _ccm )
    {
        _ccmPair.read(_cmxPair.index());
    }
    if ( _csm )
    {
        _csmPair.read(_cmxPair.index());
    }

    return _cmxPair.clusterSize();
}



/*!
 * Retrieves the sample string of the specified cluster
 * in the current gene pair. The gene pair must first
 * be set using the setPair() function.
 *
 * @param cluster_index The index of the cluster in the current gene pair.
 *
 * \return  The sample string of characters indicating the role
 * of the sample in the cluster.
 */
QString NetworkWriter::getEdgeSampleString(int cluster_index) const
{
    EDEBUG_FUNC(this,cluster_index);

    // Initialize the sample string with zeros.
    QString sampleMask(_ccm->sampleSize(), '0');

    // If cluster data exists then use it.
    if ( _ccmPair.clusterSize() > 0 )
    {
        // Write sample mask to string.
        for ( int i = 0; i < _ccm->sampleSize(); i++ )
        {
            sampleMask[i] = '0' + _ccmPair.at(cluster_index, i);
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
 * @param sample_string  The sample string as returned by the
 * getEdgeSampleString function.
 */
int NetworkWriter::getEdgeNumSamples(QString sample_string) const
{
    EDEBUG_FUNC(this,cluster_index);

    return sample_string.count("1");
}



/*!
 * If the condition-specific matrix (csm) is present then this
 * function gets the list of tests that were performed and stores
 * them in a variable local to this class for easy iteration later.
 * This function should be called in the initialize() function of
 * a child NetworkWriter implementation if the test data will be written
 * in the output file.
 */
void NetworkWriter::setTestNames()
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
 * @param cluster_index  The index of the cluster in the current gene pair.
 * @param test_index The index of the test to retrieve the value for.
 */
float NetworkWriter::getEdgeTestValue(int cluster_index, int test_index) const
{
    EDEBUG_FUNC(this,clsuster_index,test_index);

    QString fullTestName = _testNames[test_index];

    for ( int i = 0; i < _csm->getTestCount(); i++ )
    {
        QString testName = _csm->getTestName(i);

        if (fullTestName.contains("_pVal") && fullTestName.compare(testName + "_pVal") == 0)
        {
            return static_cast<float>(_csmPair.at(cluster_index, i, "pvalue"));
        }

        if (fullTestName.contains("_RSqr") && fullTestName.compare(testName + "_RSqr") == 0)
        {
            return  static_cast<float>(_csmPair.at(cluster_index, i, "r2"));
        }
    }

    return static_cast<float>(qQNaN());
}



/*!
 * Checks to make sure that the output stream is ok.
 */
void NetworkWriter::checkStatus() const
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
void FullNetworkWriter::initialize()
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
    for (int i = 0; i < _testNames.size(); ++i)
    {
        _stream << "\t" << _testNames[i];
    }

    _stream
        << "\n";
}



/*!
 * Writes a single edge corresponding to the given cluster of the current
 * gene pair to the full tab-delimited text output file. The gene pair
 * must first be set using the setPair() function.
 *
 * @param cluster_index The index of cluster in the current gene pair.
 * @param passed The set of tests that passed. (unused)
 */
void FullNetworkWriter::writeEdgeCluster(int cluster_index, QVector<QString>)
{
    EDEBUG_FUNC(this,cluster_index);

    QString source = getEdgeGene1();
    QString target = getEdgeGene2();
    float correlation = getEdgeSimilarity(cluster_index);
    QString interaction {"co"};
    int cluster_num = cluster_index + 1;
    QString sample_string = getEdgeSampleString(cluster_index);
    int num_samples = getEdgeNumSamples(sample_string);

    _stream
        << source
        << "\t" << target
        << "\t" << correlation
        << "\t" << interaction
        << "\t" << cluster_num
        << "\t" << num_samples
        << "\t" << sample_string;

    // Add each conditional test to the network file too.
    for (int i = 0; i < _testNames.size(); ++i)
    {
        float test_value = getEdgeTestValue(cluster_index, i);
        _stream << "\t" << test_value;
    }

    _stream
        << "\n";
}



/*!
 * The initialization function for the full text output. This
 * function adds the tab-delimited header to the output file.
 */
void TidyNetworkWriter::initialize()
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
 * @param cluster_index The index of cluster in the current gene pair.
 * @param passed the set of tests that passed.
 */
void TidyNetworkWriter::writeEdgeCluster(int cluster_index, QVector<QString> passed)
{
    EDEBUG_FUNC(this,cluster_index,passed);

    QString source = getEdgeGene1();
    QString target = getEdgeGene2();
    float correlation = getEdgeSimilarity(cluster_index);
    QString interaction {"co"};
    int cluster_num = cluster_index + 1;
    QString sample_string = getEdgeSampleString(cluster_index);
    int num_samples = getEdgeNumSamples(sample_string);

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
            << "\t" << num_samples
            << "\t" << sample_string
            << "\n";
    }

    // Otherwise, write out an edge for each passed test.
    else
    {
        for (int i = 0; i < passed.size(); ++i)
        {
            QString test_name = passed[i];
            int pVali = _testNames.lastIndexOf(test_name + "_pVal");
            int rSqri = _testNames.lastIndexOf(test_name + "_RSqr");
            float pVal = qQNaN();
            float rSqr = qQNaN();
            if (pVali > -1)
            {
                pVal = getEdgeTestValue(cluster_index, pVali);
            }
            if (rSqri > -1)
            {
                rSqr = getEdgeTestValue(cluster_index, rSqri);
            }
            _stream
                << source
                << "\t" << target
                << "\t" << correlation
                << "\t" << interaction
                << "\t" << cluster_num
                << "\t" << num_samples
                << "\t" << sample_string
                << "\t" << test_name
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
void MinimalNetworkWriter::initialize()
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
 * @param cluster_index The index of cluster in the current gene pair.
 * @param passed The set of tests that passed. (unused)
 */
void MinimalNetworkWriter::writeEdgeCluster(int cluster_index, QVector<QString>)
{
    EDEBUG_FUNC(this,cluster_index);

    QString source = getEdgeGene1();
    QString target = getEdgeGene2();
    float correlation = getEdgeSimilarity(cluster_index);
    int cluster_num = cluster_index + 1;

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
void GMLNetworkWriter::initialize()
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
 * @param cluster_index The index of cluster in the current gene pair.
 * @param passed The set of tests that passed. (unused)
 */
void GMLNetworkWriter::writeEdgeCluster(int cluster_index, QVector<QString>)
{
    EDEBUG_FUNC(this,cluster_index);

    QString source = getEdgeGene1();
    QString target = getEdgeGene2();
    float correlation = getEdgeSimilarity(cluster_index);
    QString interaction {"co"};
    int cluster_num = cluster_index + 1;
    QString sample_string = getEdgeSampleString(cluster_index);
    int num_samples = getEdgeNumSamples(sample_string);

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
            << "      <data key=\"cs\">" << num_samples << "</data>\n"
            << "      <data key=\"s\">" << sample_string << "</data>\n";

    // Add each conditional test to the network file too.
    for (int i = 0; i < _testNames.size(); ++i)
    {
        float test_value = getEdgeTestValue(cluster_index, i);

        _stream
            << "      <data key=\"t" << i << "\">" << test_value << "</data>\n";
    }

    // Finish out the edge.
    _stream
            << "    </edge>\n";
}



/*!
 * Terminates the GraphML output file.
 */
void GMLNetworkWriter::finish()
{
    EDEBUG_FUNC(this);

    // Write the XML header to file.
    _stream
            << "  </graph>\n"
            << "</graphml>\n";
}
