#include <ace/core/AceCore.h>
#include <ace/core/datamanager.h>
#include <ace/core/datareference.h>

#include "testcorrelationmatrix.h"
#include "correlationmatrix.h"
#include "datafactory.h"



struct Correlation
{
	GenePair::Index index;
	QVector<float> clusters;
};



void TestCorrelationMatrix::test()
{
	// create random correlation data
	int numGenes = 10;
	int maxClusters = 5;
	QVector<Correlation> testCorrelations;

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
					clusters[k] = -1.0 + 2.0 * rand() / ((2 << 31) - 1);
				}

				testCorrelations.append({ { i, j }, clusters });
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

	for ( auto& correlation : testCorrelations )
	{
		pair.clearClusters();
		pair.addCluster(correlation.clusters.size());

		for ( int k = 0; k < pair.clusterSize(); ++k )
		{
			pair.at(k, 0) = correlation.clusters.at(k);
		}

		pair.write(correlation.index);
	}

	// read and verify correlation data from file
	pair.reset();

	for ( auto& correlation : testCorrelations )
	{
		Q_ASSERT(pair.hasNext());
		pair.readNext();

		Q_ASSERT(correlation.index == pair.index());
		Q_ASSERT(pair.clusterSize() == correlation.clusters.size());

		for ( int k = 0; k < pair.clusterSize(); ++k )
		{
			Q_ASSERT(pair.at(k, 0) == correlation.clusters.at(k));
		}
	}
}
