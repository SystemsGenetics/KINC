#include <ace/core/core.h>
#include <ace/core/ace_analytic_single.h>
#include <ace/core/ace_dataobject.h>

#include "testimportcorrelationmatrix.h"
#include "../core/analyticfactory.h"
#include "../core/datafactory.h"
#include "../core/importcorrelationmatrix_input.h"



void TestImportCorrelationMatrix::test()
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
				QVector<QVector<qint8>> sampleMasks(numClusters);
				QVector<float> correlations(numClusters);

				for ( int k = 0; k < numClusters; ++k )
				{
					sampleMasks[k].resize(numSamples);

					for ( int n = 0; n < numSamples; ++n )
					{
						sampleMasks[k][n] = rand() % 2;
					}

					correlations[k] = -1.0 + 2.0 * rand() / (1 << 31);
				}

				testPairs.append({ { i, j }, sampleMasks, correlations });
			}
		}
	}

	// initialize temp files
	QString txtPath {QDir::tempPath() + "/test.txt"};
	QString ccmPath {QDir::tempPath() + "/test.ccm"};
	QString cmxPath {QDir::tempPath() + "/test.cmx"};

	QFile(txtPath).remove();
	QFile(ccmPath).remove();
	QFile(cmxPath).remove();

	// create raw text file
	QFile file(txtPath);
	QVERIFY(file.open(QIODevice::ReadWrite));

	// write correlation data to text file
	QTextStream stream(&file);

	for ( auto& testPair : testPairs )
	{
		for ( int k = 0; k < testPair.sampleMasks.size(); k++ )
		{
			int numSamples = 0;
			int numMissing = 0;
			int numPostOutliers = 0;
			int numPreOutliers = 0;
			int numThreshold = 0;
			QString sampleMask(numSamples);

			for ( int i = 0; i < numSamples; i++ )
			{
				sampleMask[i] = '0' + testPair.sampleMasks[k][i];
			}

			stream
				<< testPair.index.getX()
				<< "\t" << testPair.index.getY()
				<< "\t" << k
				<< "\t" << testPair.sampleMasks.size()
				<< "\t" << numSamples
				<< "\t" << numMissing
				<< "\t" << numPostOutliers
				<< "\t" << numPreOutliers
				<< "\t" << numThreshold
				<< "\t" << testPair.correlations[k]
				<< "\t" << sampleMask
				<< "\n";
		}
	}

	file.close();

	// create analytic manager
	auto abstractManager = Ace::Analytic::AbstractManager::makeManager(AnalyticFactory::ImportCorrelationMatrixType, 0, 1);
	auto manager = qobject_cast<Ace::Analytic::Single*>(abstractManager.release());
	manager->set(ImportCorrelationMatrix::Input::InputFile, txtPath);
	manager->set(ImportCorrelationMatrix::Input::ClusterData, ccmPath);
	manager->set(ImportCorrelationMatrix::Input::CorrelationData, cmxPath);
	manager->set(ImportCorrelationMatrix::Input::GeneSize, numGenes);
	manager->set(ImportCorrelationMatrix::Input::MaxClusterSize, maxClusters);
	manager->set(ImportCorrelationMatrix::Input::SampleSize, numSamples);
	manager->set(ImportCorrelationMatrix::Input::CorrelationName, "test");

	// run analytic
	manager->initialize();

	// TODO: wait for analytic to finish properly

	// TODO: read and verify cluster data
	// TODO: read and verify correlation data
}
