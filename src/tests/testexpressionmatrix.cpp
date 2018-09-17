#include <ace/core/core.h>
#include <ace/core/ace_dataobject.h>

#include "testexpressionmatrix.h"
#include "../core/datafactory.h"
#include "../core/expressionmatrix.h"
#include "../core/expressionmatrix_gene.h"



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
	auto transform {ExpressionMatrix::Transform::None};

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

	std::unique_ptr<Ace::DataObject> dataRef {new Ace::DataObject(path, DataFactory::ExpressionMatrixType, EMetadata(EMetadata::Object))};
	ExpressionMatrix* matrix {dataRef->data()->cast<ExpressionMatrix>()};

	// write data to file
	matrix->initialize(geneNames, sampleNames, transform);

	ExpressionMatrix::Gene gene(matrix);
	for ( int i = 0; i < matrix->geneSize(); ++i )
	{
		for ( int j = 0; j < matrix->sampleSize(); ++j )
		{
			gene[j] = testExpressions[i * numSamples + j];
		}

		gene.write(i);
	}

	matrix->finish();

	// read expression data from file
	QVector<float> expressions {matrix->dumpRawData()};

	// verify expression data
	QVERIFY(!memcmp(testExpressions.data(), expressions.data(), testExpressions.size() * sizeof(float)));
}
