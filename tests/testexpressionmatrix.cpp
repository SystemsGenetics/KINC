#include <ace/core/core.h>
#include <ace/core/ace_dataobject.h>

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
	for ( int i = 0; i < numGenes; ++i )
	{
		geneNames.append(QString::number(i));
	}

	QStringList sampleNames;
	for ( int i = 0; i < numSamples; ++i )
	{
		sampleNames.append(QString::number(i));
	}

	// create data object
	QString path {QDir::tempPath() + "/test.emx"};

	std::unique_ptr<Ace::DataObject> dataRef {new Ace::DataObject(path, DataFactory::ExpressionMatrixType, EMetadata(EMetadata::Object))};
	ExpressionMatrix* matrix {dataRef->data()->cast<ExpressionMatrix>()};

	// write data to file
	matrix->initialize(geneNames, sampleNames);

	ExpressionMatrix::Gene gene(matrix);
	for ( int i = 0; i < matrix->getGeneSize(); ++i )
	{
		for ( int j = 0; j < matrix->getSampleSize(); ++j )
		{
			gene[j] = testExpressions[i * numSamples + j];
		}

		gene.write(i);
	}

	matrix->finish();

	// read expression data from file
	std::unique_ptr<float> expressions {matrix->dumpRawData()};

	// verify expression data
	QVERIFY(!memcmp(testExpressions.data(), expressions.get(), testExpressions.size() * sizeof(float)));
}
