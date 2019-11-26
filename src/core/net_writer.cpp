#include "net_writer.h"
#include "expressionmatrix_gene.h"
#include <math.h>

NetWriter::NetWriter(QTextStream * stream, ExpressionMatrix * emx,
                     CorrelationMatrix * cmx, CCMatrix * ccm,
                     CSMatrix * csm)
{
    EDEBUG_FUNC(this);

    _stream = stream;
    _emx = emx;
    _cmx = cmx;
    _ccm = ccm;
    _csm = csm;
}

void NetWriter::initialize()
{
    EDEBUG_FUNC(this);
}

void NetWriter::writeEdgeCluster(int cluster_index)
{
    EDEBUG_FUNC(this, cluster_index);
}

void NetWriter::finish()
{
    EDEBUG_FUNC(this);
}

void NetWriter::setEdgeSampleStrings()
{
    EDEBUG_FUNC(this);

    int numSamples {0};

    for ( int k = 0; k < _cmxPair.clusterSize(); k++ )
    {

        // Initialize the sample string with zeros.
       QString sampleMask(_ccm->sampleSize(), '0');

        // If cluster data exists then use it.
        if ( _ccmPair.clusterSize() > 0 )
        {
            // Compute cluster size.
            for ( int i = 0; i < _ccm->sampleSize(); i++ )
            {
                if ( _ccmPair.at(k, i) == 1 )
                {
                    numSamples++;
                }
            }

            // Write sample mask to string.
            for ( int i = 0; i < _ccm->sampleSize(); i++ )
            {
                sampleMask[i] = '0' + _ccmPair.at(k, i);
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
                    numSamples++;
                }
            }
        }
        // Otherwise throw an error.
        else
        {
            E_MAKE_EXCEPTION(e);
            e.setTitle(QObject::tr("Invalid Input"));
            e.setDetails(QObject::tr("Expression Matrix was not provided but Cluster Matrix is missing sample data."));
            throw e;
        }

        // Save the sample string.
        _sampleStrings.append(sampleMask);
        _numSamples.append(numSamples);
    }
}

QString NetWriter::getEdgeGene1()
{
    EDEBUG_FUNC(this);

    EMetaArray geneNames {_cmx->geneNames()};
    return geneNames.at(_cmxPair.index().getX()).toString();
}

QString NetWriter::getEdgeGene2()
{
    EDEBUG_FUNC(this);

    EMetaArray geneNames {_cmx->geneNames()};
    return geneNames.at(_cmxPair.index().getY()).toString();
}

float NetWriter::getEdgeSimilarity(int cluster_index)
{
    EDEBUG_FUNC(this);

    return _cmxPair.at(cluster_index);
}

void NetWriter::setEdge(Pairwise::Index cmx_index)
{
    EDEBUG_FUNC(this);

    _index = cmx_index;
    _cmxPair.read(_index);
    _ccmPair.read(_index);
    if (_csm) {
        _csmPair.read(_index);
    }

    setEdgeSampleStrings();
}

QString NetWriter::getEdgeSampleString(int cluster_index)
{
    return _sampleStrings[cluster_index];
}

int NetWriter::getEdgeNumSamples(int cluster_index)
{
    return _numSamples[cluster_index];
}

void NetWriter::setTestNames()
{
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
float NetWriter::getEdgeTestValue(int cluster_index, int test_index)
{
    QString fullTestName = _testNames[test_index];
    for ( int i = 0; i < _csm->getTestCount(); i++ )
    {
       QString testName = _csm->getTestName(i);
       if (fullTestName.contains(testName))
       {
           if (fullTestName.contains("_pVal"))
           {
               return _csmPair.at(cluster_index, i, "pvalue");
           }
           if (fullTestName.contains("_RSqr"))
           {
               return _csmPair.at(cluster_index, i, "r2");
           }
       }
    }
    return qQNaN();
}

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

void FullNetWriter::writeEdgeCluster(int cluster_index)
{
    EDEBUG_FUNC(this);

    QString source = getEdgeGene1();
    QString target = getEdgeGene2();
    float correlation = getEdgeSimilarity(cluster_index);
    QString interaction {"co"};
    int cluster_num = cluster_index + 1;
    int num_samples = getEdgeNumSamples(cluster_index);
    QString sample_string = getEdgeSampleString(cluster_index);

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
        double test_value = getEdgeTestValue(cluster_index, i);
        *(_stream) << "\t" << test_value;
    }

    *(_stream)
        << "\n";
}

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

void MinimalNetWriter::writeEdgeCluster(int cluster_index)
{
    EDEBUG_FUNC(this);

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


void GMLNetWriter::initialize()
{
    EDEBUG_FUNC(this);



    // Write the XML header to file.
    *(_stream)
        << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n"
        << "         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
        << "         xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n"
        << "  <graph id=\"G\" edgedefault=\"undirected\">\n";

    // Write the node list to file.
    EMetaArray geneNames {_cmx->geneNames()};
    for ( int i = 0; i < _cmx->geneSize(); i++ )
    {
        QString id {geneNames.at(i).toString()};

        *(_stream) << "    <node id=\"" << id << "\"/>\n";
    }
}

void GMLNetWriter::writeEdgeCluster(int cluster_index)
{
    EDEBUG_FUNC(this);

    QString source = getEdgeGene1();
    QString target = getEdgeGene2();
    QString sample_string = getEdgeSampleString(cluster_index);

     *(_stream)
            << "    <edge"
            << "      source=\"" << source << "\""
            << "      target=\"" << target << "\""
            << "      samples=\"" << sample_string << "\""
            << "    />\n";
}


void GMLNetWriter::finish()
{
    EDEBUG_FUNC(this);


    // Write the XML header to file.
    *(_stream)
            << "  </graph>\n"
            << "</graphml>\n";

}
