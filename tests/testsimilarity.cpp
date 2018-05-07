#include <ace/core/AceCore.h>
#include <ace/core/datamanager.h>
#include <ace/core/datareference.h>

#include "testsimilarity.h"
#include "analyticfactory.h"
#include "datafactory.h"
#include "similarity.h"



void TestSimilarity::test()
{
	// create random expression data
	int numGenes = 10;
	int numSamples = 5;
	QVector<float> testExpressions(numGenes * numSamples);

	for ( int i = 0; i < testExpressions.size(); ++i )
	{
		testExpressions[i] = -10.0 + 20.0 * rand() / (1 << 31);
	}

	// create metadata
	QStringList geneNames;
	QStringList sampleNames;
	QString noSampleToken {"NA"};

	for ( int i = 0; i < numGenes; ++i )
	{
		geneNames.append(QString::number(i));
	}

	for ( int i = 0; i < numSamples; ++i )
	{
		sampleNames.append(QString::number(i));
	}

	// initialize temp files
	QString ccmPath {QDir::tempPath() + "/test.ccm"};
	QString cmxPath {QDir::tempPath() + "/test.cmx"};
	QString emxPath {QDir::tempPath() + "/test.emx"};

	QFile(ccmPath).remove();
	QFile(cmxPath).remove();
	QFile(emxPath).remove();

	// create expression matrix
	std::unique_ptr<Ace::DataReference> emxDataRef {Ace::DataManager::getInstance().open(emxPath)};
	(*emxDataRef)->clear(DataFactory::ExpressionMatrixType);

	ExpressionMatrix* emx {dynamic_cast<ExpressionMatrix*>(&((*emxDataRef)->data()))};
	emx->initialize(geneNames, sampleNames);
	emx->prepare(true);

	ExpressionMatrix::Gene gene(emx);
	for ( int i = 0; i < emx->getGeneSize(); ++i )
	{
		for ( int j = 0; j < emx->getSampleSize(); ++j )
		{
			gene[j] = testExpressions[i * numSamples + j];
		}

		gene.write(i);
	}

	emx->setTransform(ExpressionMatrix::Transform::None);

	emx->finish();
	(*emxDataRef)->writeMeta();

	// create analytic object
	EAbstractAnalyticFactory& factory {EAbstractAnalyticFactory::getInstance()};
	std::unique_ptr<EAbstractAnalytic> analytic {factory.make(AnalyticFactory::SimilarityType)};

	analytic->addDataIn(Similarity::InputData, emxPath, DataFactory::ExpressionMatrixType);
	analytic->addDataOut(Similarity::ClusterData, ccmPath, DataFactory::CCMatrixType);
	analytic->addDataOut(Similarity::CorrelationData, cmxPath, DataFactory::CorrelationMatrixType);
	analytic->setArgument(Similarity::ClusteringArg, "gmm");
	analytic->setArgument(Similarity::CorrelationArg, "pearson");

	// run analytic
	analytic->run();

	// TODO: read and verify cluster data
	// TODO: read and verify correlation data
}
