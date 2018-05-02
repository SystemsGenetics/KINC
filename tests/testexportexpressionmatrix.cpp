#include <ace/core/AceCore.h>
#include <ace/core/datamanager.h>
#include <ace/core/datareference.h>

#include "testexportexpressionmatrix.h"
#include "analyticfactory.h"
#include "datafactory.h"
#include "exportexpressionmatrix.h"



void TestExportExpressionMatrix::test()
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
	QString txtPath {QDir::tempPath() + "/test.txt"};
	QString emxPath {QDir::tempPath() + "/test.emx"};

	QFile(txtPath).remove();
	QFile(emxPath).remove();

	// create expression matrix
	std::unique_ptr<Ace::DataReference> dataRef {Ace::DataManager::getInstance().open(emxPath)};
	(*dataRef)->clear(DataFactory::ExpressionMatrixType);

	ExpressionMatrix* matrix {dynamic_cast<ExpressionMatrix*>(&((*dataRef)->data()))};
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

	matrix->setTransform(ExpressionMatrix::Transform::None);

	matrix->finish();
	(*dataRef)->writeMeta();

	// create analytic object
	EAbstractAnalyticFactory& factory {EAbstractAnalyticFactory::getInstance()};
	std::unique_ptr<EAbstractAnalytic> analytic {factory.make(AnalyticFactory::ExportExpressionMatrixType)};

	analytic->addDataIn(ExportExpressionMatrix::InputData, emxPath, DataFactory::ExpressionMatrixType);
	analytic->addFileOut(ExportExpressionMatrix::OutputFile, txtPath);
	analytic->setArgument(ExportExpressionMatrix::NoSampleToken, noSampleToken);

	// run analytic
	analytic->run();

	// read expression data from raw file
	QFile file(txtPath);
	Q_ASSERT(file.open(QIODevice::ReadOnly));

	QVector<float> expressions(numGenes * numSamples);

	QTextStream stream(&file);
	for ( int i = 0; i < numGenes + 1; ++i )
	{
		Q_ASSERT(!stream.atEnd());

		QString line = stream.readLine();
		auto words = line.splitRef(QRegExp("\\s+"),QString::SkipEmptyParts);

		if ( i == 0 )
		{
			Q_ASSERT(words.size() == numSamples);
		}
		else
		{
			Q_ASSERT(words.size() == numSamples + 1);

			for ( int j = 1; j < words.size(); ++j )
			{
				if ( words.at(j) == noSampleToken )
				{
					expressions[(i - 1) * numSamples + (j - 1)] = NAN;
				}
				else
				{
					bool ok;
					float value = words.at(j).toDouble(&ok);

					Q_ASSERT(ok);

					expressions[(i - 1) * numSamples + (j - 1)] = value;
				}
			}
		}
	}

	// verify expression data
	float error = 0;

	for ( int i = 0; i < testExpressions.size(); ++i )
	{
		error += fabs(testExpressions[i] - expressions[i]);
	}

	error /= testExpressions.size();

	Q_ASSERT(error < 1e-3);
}
