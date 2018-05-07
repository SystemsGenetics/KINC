#include <ace/core/AceCore.h>
#include <ace/core/datamanager.h>
#include <ace/core/datareference.h>

#include "testexportcorrelationmatrix.h"
#include "analyticfactory.h"
#include "datafactory.h"
#include "exportcorrelationmatrix.h"



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
	EMetadata metaSampleNames(EMetadata::Array);
	EMetadata metaCorrelationNames(EMetadata::Array);

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

	EMetadata* name {new EMetadata(EMetadata::String)};
	*(name->toString()) = "test";
	metaCorrelationNames.toArray()->append(name);

	// initialize temp files
	QString txtPath {QDir::tempPath() + "/test.txt"};
	QString ccmPath {QDir::tempPath() + "/test.ccm"};
	QString cmxPath {QDir::tempPath() + "/test.cmx"};

	QFile(txtPath).remove();
	QFile(ccmPath).remove();
	QFile(cmxPath).remove();

	// create cluster matrix
	std::unique_ptr<Ace::DataReference> ccmDataRef {Ace::DataManager::getInstance().open(ccmPath)};
	(*ccmDataRef)->clear(DataFactory::CCMatrixType);

	CCMatrix* ccm {dynamic_cast<CCMatrix*>(&((*ccmDataRef)->data()))};
	ccm->initialize(metaGeneNames, metaSampleNames);
	ccm->prepare(false);

	CCMatrix::Pair ccmPair(ccm);

	for ( auto& testPair : testPairs )
	{
		ccmPair.clearClusters();
		ccmPair.addCluster(testPair.clusters.size());

		for ( int k = 0; k < ccmPair.clusterSize(); ++k )
		{
			for ( int n = 0; n < numSamples; ++n )
			{
				ccmPair.at(k, n) = testPair.clusters.at(k).at(n);
			}
		}

		ccmPair.write(testPair.index);
	}

	ccm->finish();
	(*ccmDataRef)->writeMeta();

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
	std::unique_ptr<EAbstractAnalytic> analytic {factory.make(AnalyticFactory::ExportCorrelationMatrixType)};

	analytic->addDataIn(ExportCorrelationMatrix::ClusterData, ccmPath, DataFactory::CCMatrixType);
	analytic->addDataIn(ExportCorrelationMatrix::CorrelationData, cmxPath, DataFactory::CorrelationMatrixType);
	analytic->addFileOut(ExportCorrelationMatrix::OutputFile, txtPath);

	// run analytic
	analytic->run();

	// read and verify data from raw file
	float error = 0;
	int numClusters = 0;

	QFile file(txtPath);
	Q_ASSERT(file.open(QIODevice::ReadOnly));

	QTextStream stream(&file);

	for ( auto& testPair : testPairs )
	{
		for ( int k = 0; k < testPair.clusters.size(); ++k )
		{
			Q_ASSERT(!stream.atEnd());

			QString line = stream.readLine();
			auto words = line.splitRef(QRegExp("\\s+"), QString::SkipEmptyParts);

			Q_ASSERT(words.size() == 11);

			int geneX = words[0].toInt();
			int geneY = words[1].toInt();
			int cluster = words[2].toInt();
			int clusterSize = words[3].toInt();
			float correlation = words[9].toFloat();
			QStringRef sampleMask = words[10];

			Q_ASSERT(testPair.index.getX() == geneX);
			Q_ASSERT(testPair.index.getY() == geneY);
			Q_ASSERT(k == cluster);
			Q_ASSERT(testPair.clusters.size() == clusterSize);
			Q_ASSERT(numSamples == sampleMask.size());

			if ( testPair.clusters.size() > 1 )
			{
				for ( int i = 0; i < sampleMask.size(); ++i )
				{
					Q_ASSERT(testPair.clusters[k][i] == sampleMask[i].digitValue());
				}
			}

			error += fabs(testPair.correlations[k] - correlation);
			numClusters++;
		}
	}

	error /= numClusters;

	Q_ASSERT(error < 1e-3);
}
