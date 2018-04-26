#include <ace/core/AceCore.h>
#include <ace/core/datamanager.h>
#include <ace/core/datareference.h>

#include "testexpressionmatrix.h"
#include "datafactory.h"
#include "expressionmatrix.h"



void TestExpressionMatrix::test()
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

	for ( int i = 0; i < numGenes; ++i )
	{
		geneNames.append(QString::number(i));
	}

	for ( int i = 0; i < numSamples; ++i )
	{
		sampleNames.append(QString::number(i));
	}

	// create data object
	QString path {QDir::tempPath() + "/test.emx"};

	std::unique_ptr<Ace::DataReference> dataRef {Ace::DataManager::getInstance().open(path)};
   (*dataRef)->clear(DataFactory::ExpressionMatrixType);

	ExpressionMatrix* matrix {dynamic_cast<ExpressionMatrix*>(&((*dataRef)->data()))};

	// write data to file
	matrix->initialize(geneNames, sampleNames);
	matrix->prepare(true);

	ExpressionMatrix::Gene gene(matrix);

	for ( int i = 0; i < matrix->getGeneSize(); ++i )
	{
		for ( int j = 0; j < matrix->getSampleSize(); ++j )
		{
			gene[j] = testExpressions[i * numSamples + j];
		}

		gene.write(i);
	}

	// read expression data from file
	std::unique_ptr<float> expressions {matrix->dumpRawData()};

	// verify expression data
	Q_ASSERT(!memcmp(testExpressions.data(), expressions.get(), testExpressions.size() * sizeof(float)));
}
