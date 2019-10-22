#include <ace/core/core.h>
#include <ace/core/ace_analytic_single.h>
#include <ace/core/ace_dataobject.h>

#include "testexportcorrelationmatrix.h"
#include "../core/analyticfactory.h"
#include "../core/datafactory.h"
#include "../core/exportcorrelationmatrix_input.h"
#include "../core/ccmatrix_pair.h"
#include "../core/correlationmatrix_pair.h"



void TestExportCorrelationMatrix::test()
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

					correlations[k] = -1.0f + 2.0f * rand() / (1 << 31);
				}

				testPairs.append({ { i, j }, sampleMasks, correlations });
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

	QString correlationName("test");

	// initialize temp files
	QString txtPath {QDir::tempPath() + "/test.txt"};
	QString ccmPath {QDir::tempPath() + "/test.ccm"};
	QString cmxPath {QDir::tempPath() + "/test.cmx"};

	QFile(txtPath).remove();
	QFile(ccmPath).remove();
	QFile(cmxPath).remove();

	// create cluster matrix
	std::unique_ptr<Ace::DataObject> ccmDataRef {new Ace::DataObject(ccmPath)};
	CCMatrix* ccm {ccmDataRef->data()->cast<CCMatrix>()};

	ccm->initialize(metaGeneNames, maxClusters, metaSampleNames);

	CCMatrix::Pair ccmPair(ccm);
	for ( auto& testPair : testPairs )
	{
		ccmPair.clearClusters();
		ccmPair.addCluster(testPair.sampleMasks.size());

		for ( int k = 0; k < ccmPair.clusterSize(); ++k )
		{
			for ( int n = 0; n < numSamples; ++n )
			{
				ccmPair.at(k, n) = testPair.sampleMasks.at(k).at(n);
			}
		}

		ccmPair.write(testPair.index);
	}

	ccmDataRef->data()->finish();
	ccmDataRef->finalize();

	// create correlation matrix
	std::unique_ptr<Ace::DataObject> cmxDataRef {new Ace::DataObject(cmxPath)};
	CorrelationMatrix* cmx {cmxDataRef->data()->cast<CorrelationMatrix>()};

	cmx->initialize(metaGeneNames, maxClusters, correlationName);

	CorrelationMatrix::Pair cmxPair(cmx);
	for ( auto& testPair : testPairs )
	{
		cmxPair.clearClusters();
		cmxPair.addCluster(testPair.correlations.size());

		for ( int k = 0; k < cmxPair.clusterSize(); ++k )
		{
			cmxPair.at(k) = testPair.correlations.at(k);
		}

		cmxPair.write(testPair.index);
	}

	cmxDataRef->data()->finish();
	cmxDataRef->finalize();

	// create analytic manager
	auto abstractManager = Ace::Analytic::AbstractManager::makeManager(AnalyticFactory::ExportCorrelationMatrixType, 0, 1);
	auto manager = qobject_cast<Ace::Analytic::Single*>(abstractManager.release());
	manager->set(ExportCorrelationMatrix::Input::ClusterData, ccmPath);
	manager->set(ExportCorrelationMatrix::Input::CorrelationData, cmxPath);
	manager->set(ExportCorrelationMatrix::Input::OutputFile, txtPath);

	// run analytic
	manager->initialize();

	// TODO: wait for analytic to finish properly

	// read and verify data from raw file
	float error = 0;
	int numClusters = 0;

	QFile file(txtPath);
	QVERIFY(file.open(QIODevice::ReadOnly));

	QTextStream stream(&file);
	for ( auto& testPair : testPairs )
	{
		for ( int k = 0; k < testPair.sampleMasks.size(); ++k )
		{
			QVERIFY(!stream.atEnd());

			QString line = stream.readLine();
			auto words = line.splitRef(QRegExp("\\s+"), QString::SkipEmptyParts);

			QCOMPARE(words.size(), 7);

			int geneX = words[0].toInt();
			int geneY = words[1].toInt();
			int cluster = words[2].toInt();
			int clusterSize = words[3].toInt();
			float correlation = words[5].toFloat();
			QStringRef sampleMask = words[6];

			QCOMPARE(geneX, testPair.index.getX());
			QCOMPARE(geneY, testPair.index.getY());
			QCOMPARE(cluster, k);
			QCOMPARE(clusterSize, testPair.sampleMasks.size());
			QCOMPARE(sampleMask.size(), numSamples);

			if ( testPair.sampleMasks.size() > 1 )
			{
				for ( int i = 0; i < sampleMask.size(); ++i )
				{
					QCOMPARE((qint8) sampleMask[i].digitValue(), testPair.sampleMasks[k][i]);
				}
			}

			error += fabsf(testPair.correlations[k] - correlation);
			numClusters++;
		}
	}

	error /= numClusters;

	QVERIFY(error < 1e-3f);
}
