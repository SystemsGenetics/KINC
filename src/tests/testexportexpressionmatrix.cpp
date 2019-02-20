#include <ace/core/core.h>
#include <ace/core/ace_analytic_single.h>
#include <ace/core/ace_dataobject.h>

#include "testexportexpressionmatrix.h"
#include "../core/analyticfactory.h"
#include "../core/datafactory.h"
#include "../core/exportexpressionmatrix_input.h"
#include "../core/expressionmatrix_gene.h"



void TestExportExpressionMatrix::test()
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
	QString nanToken {"NA"};

	for ( int i = 0; i < numGenes; ++i )
	{
		geneNames.append(QString::number(i));
	}

	for ( int i = 0; i < numSamples; ++i )
	{
		sampleNames.append(QString::number(i));
	}

	// initialize temp files
	QString txtPath {QDir::tempPath() + "/test.txt"};
	QString emxPath {QDir::tempPath() + "/test.emx"};

	QFile(txtPath).remove();
	QFile(emxPath).remove();

	// create expression matrix
	std::unique_ptr<Ace::DataObject> dataRef {new Ace::DataObject(emxPath)};
	ExpressionMatrix* matrix {dataRef->data()->cast<ExpressionMatrix>()};

	matrix->initialize(geneNames, sampleNames);

	ExpressionMatrix::Gene gene(matrix);
	for ( int i = 0; i < matrix->geneSize(); ++i )
	{
		for ( int j = 0; j < matrix->sampleSize(); ++j )
		{
			gene[j] = testExpressions[i * numSamples + j];
		}

		gene.write(i);
	}

	dataRef->data()->finish();
	dataRef->finalize();

	// create analytic manager
	auto abstractManager = Ace::Analytic::AbstractManager::makeManager(AnalyticFactory::ExportExpressionMatrixType, 0, 1);
	auto manager = qobject_cast<Ace::Analytic::Single*>(abstractManager.release());
	manager->set(ExportExpressionMatrix::Input::InputData, emxPath);
	manager->set(ExportExpressionMatrix::Input::OutputFile, txtPath);
	manager->set(ExportExpressionMatrix::Input::NANToken, nanToken);

	// run analytic
	manager->initialize();

	// TODO: wait for analytic to finish properly

	// read expression data from raw file
	QFile file(txtPath);
	QVERIFY(file.open(QIODevice::ReadOnly));

	QVector<float> expressions(numGenes * numSamples);

	QTextStream stream(&file);
	for ( int i = 0; i < numGenes + 1; ++i )
	{
		QVERIFY(!stream.atEnd());

		QString line = stream.readLine();
		auto words = line.splitRef(QRegExp("\\s+"),QString::SkipEmptyParts);

		if ( i == 0 )
		{
			QCOMPARE(words.size(), numSamples);
		}
		else
		{
			QCOMPARE(words.size(), numSamples + 1);

			for ( int j = 1; j < words.size(); ++j )
			{
				if ( words.at(j) == nanToken )
				{
					expressions[(i - 1) * numSamples + (j - 1)] = NAN;
				}
				else
				{
					bool ok;
					float value = words.at(j).toFloat(&ok);

					QVERIFY(ok);

					expressions[(i - 1) * numSamples + (j - 1)] = value;
				}
			}
		}
	}

	// verify expression data
	float error = 0;

	for ( int i = 0; i < testExpressions.size(); ++i )
	{
		error += fabsf(testExpressions[i] - expressions[i]);
	}

	error /= testExpressions.size();

	QVERIFY(error < 1e-3f);
}
