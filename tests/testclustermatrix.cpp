#include <ace/core/AceCore.h>
#include <ace/core/datamanager.h>
#include <ace/core/datareference.h>

#include "testclustermatrix.h"
#include "ccmatrix.h"
#include "datafactory.h"



void TestClusterMatrix::test()
{
	// create random cluster data
	int numGenes = 10;
	int numSamples = 5;
	int maxClusters = 5;
	QVector<Pair> testPairs;

	for ( int i = 0; i < numGenes; ++i )
	{
		for ( int j = 0; j < i; ++j )
		{
			int numClusters = rand() % (maxClusters + 1);

			if ( numClusters > 0 )
			{
				QVector<QVector<qint8>> clusters(numClusters);

				for ( int k = 0; k < numClusters; ++k )
				{
					clusters[k].resize(numSamples);

					for ( int n = 0; n < numSamples; ++n )
					{
						clusters[k][n] = rand() % 16;
					}
				}

				testPairs.append({ { i, j }, clusters });
			}
		}
	}

	// create metadata
	EMetadata metaGeneNames(EMetadata::Array);
	EMetadata metaSampleNames(EMetadata::Array);

	for ( int i = 0; i < numGenes; ++i )
	{
		EMetadata* name {new EMetadata(EMetadata::String)};
		*(name->toString()) = QString::number(i);
		metaGeneNames.toArray()->append(name);
	}

	for ( int i = 0; i < numSamples; ++i )
	{
		EMetadata* name {new EMetadata(EMetadata::String)};
		*(name->toString()) = QString::number(i);
		metaSampleNames.toArray()->append(name);
	}

	// create data object
	QString path {QDir::tempPath() + "/test.ccm"};

	std::unique_ptr<Ace::DataReference> dataRef {Ace::DataManager::getInstance().open(path)};
	(*dataRef)->clear(DataFactory::CCMatrixType);

	CCMatrix* matrix {dynamic_cast<CCMatrix*>(&((*dataRef)->data()))};

	// write data to file
	matrix->initialize(metaGeneNames, metaSampleNames);
	matrix->prepare(false);

	CCMatrix::Pair pair(matrix);

	for ( auto& testPair : testPairs )
	{
		pair.clearClusters();
		pair.addCluster(testPair.clusters.size());

		for ( int k = 0; k < pair.clusterSize(); ++k )
		{
			for ( int n = 0; n < numSamples; ++n )
			{
				pair.at(k, n) = testPair.clusters.at(k).at(n);
			}
		}

		pair.write(testPair.index);
	}

	matrix->finish();
	(*dataRef)->writeMeta();

	// read and verify correlation data from file
	pair.reset();

	for ( auto& testPair : testPairs )
	{
		Q_ASSERT(pair.hasNext());
		pair.readNext();

		Q_ASSERT(testPair.index == pair.index());
		Q_ASSERT(pair.clusterSize() == testPair.clusters.size());

		for ( int k = 0; k < pair.clusterSize(); ++k )
		{
			for ( int n = 0; n < numSamples; ++n )
			{
				Q_ASSERT(pair.at(k, n) == testPair.clusters.at(k).at(n));
			}
		}
	}
}
