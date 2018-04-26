#include <ace/core/AceCore.h>
#include <ace/core/datamanager.h>
#include <ace/core/datareference.h>

#include "testcorrelationmatrix.h"
#include "correlationmatrix.h"
#include "datafactory.h"



struct CorrelationPair
{
	GenePair::Index index;
	QVector<float> clusters;
};



void TestCorrelationMatrix::test()
{
	// create random correlation data
	int numGenes = 10;
	int maxClusters = 5;
	QVector<CorrelationPair> testPairs;

	for ( int i = 0; i < numGenes; ++i )
	{
		for ( int j = 0; j < i; ++j )
		{
			int numClusters = rand() % (maxClusters + 1);

			if ( numClusters > 0 )
			{
				QVector<float> clusters(numClusters);

				for ( int k = 0; k < numClusters; ++k )
				{
					clusters[k] = -1.0 + 2.0 * rand() / (1 << 31);
				}

				testPairs.append({ { i, j }, clusters });
			}
		}
	}

	// create metadata
   EMetadata metaGeneNames(EMetadata::Array);
   EMetadata metaCorrelationNames(EMetadata::Array);

   for ( int i = 0; i < numGenes; ++i )
   {
      EMetadata* name {new EMetadata(EMetadata::String)};
      *(name->toString()) = QString::number(i);
      metaGeneNames.toArray()->append(name);
   }

   EMetadata* name {new EMetadata(EMetadata::String)};
   *(name->toString()) = "test";
   metaCorrelationNames.toArray()->append(name);

	// create data object
	QString path {QDir::tempPath() + "/test.cmx"};

	std::unique_ptr<Ace::DataReference> dataRef {Ace::DataManager::getInstance().open(path)};
   (*dataRef)->clear(DataFactory::CorrelationMatrixType);

	CorrelationMatrix* matrix {dynamic_cast<CorrelationMatrix*>(&((*dataRef)->data()))};

	// write data to file
	matrix->initialize(metaGeneNames, metaCorrelationNames);
	matrix->prepare(false);

	CorrelationMatrix::Pair pair(matrix);

	for ( auto& testPair : testPairs )
	{
		pair.clearClusters();
		pair.addCluster(testPair.clusters.size());

		for ( int k = 0; k < pair.clusterSize(); ++k )
		{
			pair.at(k, 0) = testPair.clusters.at(k);
		}

		pair.write(testPair.index);
	}

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
			Q_ASSERT(pair.at(k, 0) == testPair.clusters.at(k));
		}
	}
}
