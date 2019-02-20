#include <ace/core/core.h>
#include <ace/core/ace_analytic_single.h>
#include <ace/core/ace_dataobject.h>

#include "testrmt.h"
#include "../core/analyticfactory.h"
#include "../core/datafactory.h"
#include "../core/rmt_input.h"
#include "../core/correlationmatrix.h"
#include "../core/correlationmatrix_pair.h"



void TestRMT::test()
{
	// create random correlation data
	int numGenes = 10;
	int maxClusters = 5;
	QVector<Pair> testPairs;

	for ( int i = 0; i < numGenes; ++i )
	{
		for ( int j = 0; j < i; ++j )
		{
			int numClusters = rand() % (maxClusters + 1);

			if ( numClusters > 0 )
			{
				QVector<float> correlations(numClusters);

				for ( int k = 0; k < numClusters; ++k )
				{
					correlations[k] = -1.0f + 2.0f * rand() / (1 << 31);
				}

				testPairs.append({ { i, j }, correlations });
			}
		}
	}

	// create metadata
	EMetaArray metaGeneNames;
	for ( int i = 0; i < numGenes; ++i )
	{
		metaGeneNames.append(QString::number(i));
	}

	EMetaArray metaCorrelationNames;
	metaCorrelationNames.append(QString("test"));

	// initialize temp files
	QString cmxPath {QDir::tempPath() + "/test.cmx"};
	QString logPath {QDir::tempPath() + "/test.log"};

	QFile(cmxPath).remove();
	QFile(logPath).remove();

	// create correlation matrix
	std::unique_ptr<Ace::DataObject> cmxDataRef {new Ace::DataObject(cmxPath)};
	CorrelationMatrix* cmx {cmxDataRef->data()->cast<CorrelationMatrix>()};

	cmx->initialize(metaGeneNames, maxClusters, metaCorrelationNames);

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

	cmxDataRef->data()->finish();
	cmxDataRef->finalize();

	// create analytic manager
	auto abstractManager = Ace::Analytic::AbstractManager::makeManager(AnalyticFactory::RMTType, 0, 1);
	auto manager = qobject_cast<Ace::Analytic::Single*>(abstractManager.release());
	manager->set(RMT::Input::InputData, cmxPath);
	manager->set(RMT::Input::LogFile, logPath);

	// run analytic
	manager->initialize();

	// TODO: wait for analytic to finish properly
}
