#include <ace/core/core.h>
#include <ace/core/ace_dataobject.h>

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
				QVector<QVector<qint8>> sampleMasks(numClusters);

				for ( int k = 0; k < numClusters; ++k )
				{
					sampleMasks[k].resize(numSamples);

					for ( int n = 0; n < numSamples; ++n )
					{
						sampleMasks[k][n] = rand() % 16;
					}
				}

				testPairs.append({ { i, j }, sampleMasks });
			}
		}
	}

	// create metadata
	EMetaArray metaGeneNames;
	for ( int i = 0; i < numGenes; ++i )
	{
		metaGeneNames.append(QString::number(i));
	}

	EMetaArray metaSampleNames;
	for ( int i = 0; i < numSamples; ++i )
	{
		metaSampleNames.append(QString::number(i));
	}

	// create data object
	QString path {QDir::tempPath() + "/test.ccm"};

	std::unique_ptr<Ace::DataObject> dataRef {new Ace::DataObject(path, DataFactory::CCMatrixType, EMetadata(EMetadata::Object))};
	CCMatrix* matrix {qobject_cast<CCMatrix*>(dataRef->data())};

	// write data to file
	matrix->initialize(metaGeneNames, maxClusters, metaSampleNames);

	CCMatrix::Pair pair(matrix);

	for ( auto& testPair : testPairs )
	{
		pair.clearClusters();
		pair.addCluster(testPair.sampleMasks.size());

		for ( int k = 0; k < pair.clusterSize(); ++k )
		{
			for ( int n = 0; n < numSamples; ++n )
			{
				pair.at(k, n) = testPair.sampleMasks.at(k).at(n);
			}
		}

		pair.write(testPair.index);
	}

	matrix->finish();

	// read and verify correlation data from file
	pair.reset();

	for ( auto& testPair : testPairs )
	{
		QVERIFY(pair.hasNext());
		pair.readNext();

		QCOMPARE(pair.index(), testPair.index);
		QCOMPARE(pair.clusterSize(), testPair.sampleMasks.size());

		for ( int k = 0; k < pair.clusterSize(); ++k )
		{
			for ( int n = 0; n < numSamples; ++n )
			{
				QCOMPARE(pair.at(k, n), testPair.sampleMasks.at(k).at(n));
			}
		}
	}
}
