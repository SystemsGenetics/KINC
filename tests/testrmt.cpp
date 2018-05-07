#include <ace/core/AceCore.h>
#include <ace/core/datamanager.h>
#include <ace/core/datareference.h>

#include "testrmt.h"
#include "analyticfactory.h"
#include "datafactory.h"
#include "rmt.h"
#include "correlationmatrix.h"



void TestRMT::test()
{
	// create random correlation data
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
				QVector<float> correlations(numClusters);

				for ( int k = 0; k < numClusters; ++k )
				{
					clusters[k].resize(numSamples);

					for ( int n = 0; n < numSamples; ++n )
					{
						clusters[k][n] = rand() % 2;
					}

					correlations[k] = -1.0 + 2.0 * rand() / (1 << 31);
				}

				testPairs.append({ { i, j }, clusters, correlations });
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

	// initialize temp files
	QString cmxPath {QDir::tempPath() + "/test.cmx"};
	QString logPath {QDir::tempPath() + "/test.log"};

	QFile(cmxPath).remove();
	QFile(logPath).remove();

	// create correlation matrix
	std::unique_ptr<Ace::DataReference> cmxDataRef {Ace::DataManager::getInstance().open(cmxPath)};
	(*cmxDataRef)->clear(DataFactory::CorrelationMatrixType);

	CorrelationMatrix* cmx {dynamic_cast<CorrelationMatrix*>(&((*cmxDataRef)->data()))};
	cmx->initialize(metaGeneNames, metaCorrelationNames);
	cmx->prepare(false);

	CorrelationMatrix::Pair cmxPair(cmx);

	for ( auto& testPair : testPairs )
	{
		cmxPair.clearClusters();
		cmxPair.addCluster(testPair.correlations.size());

		for ( int k = 0; k < cmxPair.clusterSize(); ++k )
		{
			cmxPair.at(k, 0) = testPair.correlations.at(k);
		}

		cmxPair.write(testPair.index);
	}

	cmx->finish();
	(*cmxDataRef)->writeMeta();

	// create analytic object
	EAbstractAnalyticFactory& factory {EAbstractAnalyticFactory::getInstance()};
	std::unique_ptr<EAbstractAnalytic> analytic {factory.make(AnalyticFactory::RMTType)};

	analytic->addDataIn(RMT::InputData, cmxPath, DataFactory::CorrelationMatrixType);
	analytic->addFileOut(RMT::LogFile, logPath);

	// run analytic
	analytic->run();
}
