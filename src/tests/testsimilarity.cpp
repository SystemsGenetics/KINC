#include <ace/core/core.h>
#include <ace/core/ace_analytic_single.h>
#include <ace/core/ace_dataobject.h>

#include "testsimilarity.h"
#include "../core/analyticfactory.h"
#include "../core/datafactory.h"
#include "../core/similarity_input.h"
#include "../core/expressionmatrix_gene.h"



void TestSimilarity::test()
{
	// create random expression data
	int numGenes = 10;
	int numSamples = 5;
	QVector<float> testExpressions(numGenes * numSamples);

	for ( int i = 0; i < testExpressions.size(); ++i )
	{
		testExpressions[i] = -10.0f + 20.0f * rand() / (1 << 31);
	}

	// create metadata
	QStringList geneNames;
	QStringList sampleNames;

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
	std::unique_ptr<Ace::DataObject> emxDataRef {new Ace::DataObject(emxPath)};
	ExpressionMatrix* emx {emxDataRef->data()->cast<ExpressionMatrix>()};

	emx->initialize(geneNames, sampleNames);

	ExpressionMatrix::Gene gene(emx);
	for ( int i = 0; i < emx->geneSize(); ++i )
	{
		for ( int j = 0; j < emx->sampleSize(); ++j )
		{
			gene[j] = testExpressions[i * numSamples + j];
		}

		gene.write(i);
	}

	emxDataRef->data()->finish();
	emxDataRef->finalize();

	// create analytic manager
	auto abstractManager = Ace::Analytic::AbstractManager::makeManager(AnalyticFactory::SimilarityType, 0, 1);
	auto manager = qobject_cast<Ace::Analytic::Single*>(abstractManager.release());
	manager->set(Similarity::Input::InputData, emxPath);
	manager->set(Similarity::Input::ClusterData, ccmPath);
	manager->set(Similarity::Input::CorrelationData, cmxPath);
	manager->set(Similarity::Input::ClusteringType, "gmm");
	manager->set(Similarity::Input::CorrelationType, "pearson");

	// run analytic
	manager->initialize();

	// TODO: wait for analytic to finish properly

	// TODO: read and verify cluster data
	// TODO: read and verify correlation data
}
