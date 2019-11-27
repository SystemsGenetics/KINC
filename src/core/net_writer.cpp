#include "net_writer.h"
#include "expressionmatrix_gene.h"
#include <math.h>

/*!
 * The constructor for the NetWriter class.
 *
 * \param stream  A stream object pointer to which the output
 *        will be writtine.
 * \param emx  The expression matrix.
 * \param cmx  The correlation matrix.
 * \param ccm  The cluster composition matrix.
 * \param csm  The condition-specific matrix.
 */
NetWriter::NetWriter(QTextStream * stream, ExpressionMatrix * emx,
                     CorrelationMatrix * cmx, CCMatrix * ccm,
                     CSMatrix * csm)
{
    EDEBUG_FUNC(this, stream, emx, cmx, ccm, csm);

    _stream = stream;
    _emx = emx;

    _cmx = cmx;
    _cmxPair = CorrelationMatrix::Pair(_cmx);

    _ccm = ccm;
    _ccmPair = CCMatrix::Pair(_ccm);

    _csm = csm;
    _csmPair = CSMatrix::Pair(_csm);
}





/*!
 * Returns the name of the first gene in a gene pair. The gene pair must first
 * be set using the setPair() function.
 */
QString NetWriter::getEdgeGene1()
{
    EDEBUG_FUNC(this);

    EMetaArray geneNames {_cmx->geneNames()};
    return geneNames.at(_cmxPair.index().getX()).toString();
}





/*!
 * Returns the name of the second gene in a gene pair. The gene pair must first
 * be set using the setPair() function.
 */
QString NetWriter::getEdgeGene2()
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
 * \param cluster_index  The index of the cluster in the gene pair.
 */
float NetWriter::getEdgeSimilarity(int cluster_index)
{
    EDEBUG_FUNC(this, cluster_index);

    return _cmxPair.at(cluster_index);
}





/*!
 * Sets the current gene pair.
 *
 * \param cmx_index  The index of the pair from the correlation
 * matrix.  The correlation matrix is used to ensure all other
 * matricies (ccm, csm) are all on the same pair.
 */
void NetWriter::setPair(Pairwise::Index cmx_index)
{
    EDEBUG_FUNC(this,  cmx_index);

    _index = cmx_index;
    _cmxPair.read(_index);
    _ccmPair.read(_index);
    if (_csm) {
        _csmPair.read(_index);
    }
}





/*!
 * Retrieves the sample string of the specified cluster
 * in the current gene pair. The gene pair must first
 * be set using the setPair() function.
 *
 * \param cluster_index The index of the cluster in the current gene pair.
 *
 * \return  The sample string of characters indicating the role
 * of the sample in the cluster.
 */
QString NetWriter::getEdgeSampleString(int cluster_index)
{
    EDEBUG_FUNC(this, cluster_index);

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
 * \param sample_string  The sample string as returned by the
 * getEdgeSampleString function.
 */
int NetWriter::getEdgeNumSamples(QString sample_string)
{
    EDEBUG_FUNC(this, cluster_index);
    return sample_string.count("1");
}





/*!
 * If the condition-specific matrix (csm) is present then this
 * function gets the list of tests that were performed and stores
 * them in a variable local to this class for easy iteration later.
 * This function should be called in the initialize() function of
 * a child NetWriter implementation if the test data will be written
 * in the output file.
 */
void NetWriter::setTestNames()
{
    EDEBUG_FUNC(this);

    if (!_csm) {
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
 * \param cluster_index  The index of the cluster in the current gene pair.
 * \param test_index The index of the test to retrieve the value for.
 */
float NetWriter::getEdgeTestValue(int cluster_index, int test_index)
{
    EDEBUG_FUNC(this, clsuster_index, test_index);

    QString fullTestName = _testNames[test_index];
    for ( int i = 0; i < _csm->getTestCount(); i++ )
    {
       QString testName = _csm->getTestName(i);
       if (fullTestName.contains(testName))
       {
           if (fullTestName.contains("_pVal"))
           {
               return static_cast<float>(_csmPair.at(cluster_index, i, "pvalue"));
           }
           if (fullTestName.contains("_RSqr"))
           {
               return  static_cast<float>(_csmPair.at(cluster_index, i, "r2"));
           }
       }
    }
    return static_cast<float>(qQNaN());
}





/*!
 * The initialization function for the full text output. This
 * function adds the tab-delimited header to the output file.
 */
void FullNetWriter::initialize()
{
    EDEBUG_FUNC(this);

    setTestNames();

    // Start by writing the tab-delimited header.
    *(_stream)
        << "Source"
        << "\t" << "Target"
        << "\t" << "Similarity_Score"
        << "\t" << "Interaction"
        << "\t" << "Cluster_Index"
        << "\t" << "Cluster_Size"
        << "\t" << "Samples";

    // Add each conditional test to the network file too.
    for (int i = 0; i < _testNames.size(); ++i) {
        *(_stream) << "\t" << _testNames[i];
    }

    *(_stream)
        << "\n";
}





/*!
 * Writes a single edge corresponding to the given cluster of the current
 * gene pair to the full tab-delimited text output file. The gene pair
 * must first be set using the setPair() function.
 *
 * \param cluster_index The index of cluster in the current gene pair.
 */
void FullNetWriter::writeEdgeCluster(int cluster_index)
{
    EDEBUG_FUNC(this, cluster_index);

    QString source = getEdgeGene1();
    QString target = getEdgeGene2();
    float correlation = getEdgeSimilarity(cluster_index);
    QString interaction {"co"};
    int cluster_num = cluster_index + 1;
    QString sample_string = getEdgeSampleString(cluster_index);
    int num_samples = getEdgeNumSamples(sample_string);

     *(_stream)
            << source
            << "\t" << target
            << "\t" << correlation
            << "\t" << interaction
            << "\t" << cluster_num
            << "\t" << num_samples
            << "\t" << sample_string;

    // Add each conditional test to the network file too.
    for (int i = 0; i < _testNames.size(); ++i) {
        float test_value = getEdgeTestValue(cluster_index, i);
        *(_stream) << "\t" << test_value;
    }

    *(_stream)
        << "\n";
}





/*!
 * The initialization function for the minimial text output. This
 * function adds the tab-delimited header to the output file.
 */
void MinimalNetWriter::initialize()
{
    EDEBUG_FUNC(this);

    // Start by writing the tab-delimited header.
    *(_stream)
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
 * \param cluster_index The index of cluster in the current gene pair.
 */
void MinimalNetWriter::writeEdgeCluster(int cluster_index)
{
    EDEBUG_FUNC(this, cluster_index);

    QString source = getEdgeGene1();
    QString target = getEdgeGene2();
    float correlation = getEdgeSimilarity(cluster_index);
    int cluster_num = cluster_index + 1;

     *(_stream)
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
void GMLNetWriter::initialize()
{
    EDEBUG_FUNC(this);

    setTestNames();

    // Write the XML header to file.
    *(_stream)
        << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n"
        << "         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
        << "         xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n";

    // Add the attributes definitions
    *(_stream)
        << "  <key id=\"sc\" for=\"edge\" attr.name=\"Similarity_Score\" attr.type=\"string\"/>\n"
        << "  <key id=\"int\" for=\"edge\" attr.name=\"Interaction\" attr.type=\"string\"/>\n"
        << "  <key id=\"ci\" for=\"edge\" attr.name=\"Cluster_Index\" attr.type=\"int\"/>\n"
        << "  <key id=\"cs\" for=\"edge\" attr.name=\"Cluster_Size\" attr.type=\"int\"/>\n"
        << "  <key id=\"s\" for=\"edge\" attr.name=\"Samples\" attr.type=\"string\"/>\n";

    // Add attributes for the tests.
    for (int i = 0; i < _testNames.size(); ++i) {
        *(_stream)
          << "  <key id=\"t" << i << "\" for=\"node\" attr.name=\"" << _testNames[i] << "\" attr.type=\"double\"/>\n";
    }

    // Start the graph.
    *(_stream)
        << "  <graph id=\"G01\" edgedefault=\"undirected\">\n";
}




/*!
 * Writes a single edge corresponding to the given cluster of the current
 * gene pair to the GraphML output file. The gene pair
 * must first be set using the setPair() function.
 *
 * \param cluster_index The index of cluster in the current gene pair.
 */
void GMLNetWriter::writeEdgeCluster(int cluster_index)
{
    EDEBUG_FUNC(this, cluster_index);

    QString source = getEdgeGene1();
    QString target = getEdgeGene2();
    float correlation = getEdgeSimilarity(cluster_index);
    QString interaction {"co"};
    int cluster_num = cluster_index + 1;
    QString sample_string = getEdgeSampleString(cluster_index);
    int num_samples = getEdgeNumSamples(sample_string);

    // If we've never seen this node before then add it to
    // our list and to the output file.
    if (!_nodes.contains(source)) {
        _nodes.insert(source, true);
        *(_stream)
            << "    <node id=\"" << source << "\"/>\n";
    }
    if (!_nodes.contains(target)) {
        _nodes.insert(target, true);
        *(_stream)
            << "    <node id=\"" << target << "\"/>\n";
    }

    // Start the edge and add in the attributes.
     *(_stream)
            << "    <edge source=\"" << source << "\" target=\"" << target << "\">\n"
            << "      <data key=\"sc\">" << correlation << "</data>\n"
            << "      <data key=\"int\">co</data>\n"
            << "      <data key=\"ci\">" << cluster_num << "</data>\n"
            << "      <data key=\"cs\">" << num_samples << "</data>\n"
            << "      <data key=\"s\">" << sample_string << "</data>\n";

    // Add each conditional test to the network file too.
    for (int i = 0; i < _testNames.size(); ++i) {
        float test_value = getEdgeTestValue(cluster_index, i);
        *(_stream)
                << "      <data key=\"t" << i << "\">" << test_value << "</data>\n";
    }

    // Finish out the edge.
    *(_stream)
            << "    </edge>\n";

}




/*!
 * Terminates the GraphML output file.
 */
void GMLNetWriter::finish()
{
    EDEBUG_FUNC(this);

    // Write the XML header to file.
    *(_stream)
            << "  </graph>\n"
            << "</graphml>\n";

}
