#include <ace/core/AceCore.h>
#include <ace/core/datamanager.h>
#include <ace/core/datareference.h>

#include "testimportexpressionmatrix.h"
#include "analyticfactory.h"
#include "datafactory.h"
#include "importexpressionmatrix.h"



void TestImportExpressionMatrix::test()
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

	// create raw text file
	QFile file(txtPath);
	Q_ASSERT(file.open(QIODevice::ReadWrite));

	QTextStream stream(&file);

	for ( int i = 0; i < numSamples; i++ )
	{
		stream << sampleNames.at(i) << "\t";
	}
	stream << "\n";

	for ( int i = 0; i < numGenes; i++ )
	{
		stream << geneNames.at(i);

		for ( int j = 0; j < numSamples; j++ )
		{
			float value = testExpressions[i * numSamples + j];

			if ( std::isnan(value) )
			{
				stream << "\t" << noSampleToken;
			}
			else
			{
				stream << "\t" << value;
			}
		}

		stream << "\n";
	}

	file.close();

	// create analytic object
	EAbstractAnalyticFactory& factory {EAbstractAnalyticFactory::getInstance()};
	std::unique_ptr<EAbstractAnalytic> analytic {factory.make(AnalyticFactory::ImportExpressionMatrixType)};

	analytic->addFileIn(ImportExpressionMatrix::InputFile, txtPath);
	analytic->addDataOut(ImportExpressionMatrix::OutputData, emxPath, DataFactory::ExpressionMatrixType);
	analytic->setArgument(ImportExpressionMatrix::NoSampleToken, noSampleToken);

	// run analytic
	analytic->run();

	// read expression data from file
	std::unique_ptr<Ace::DataReference> dataRef {Ace::DataManager::getInstance().open(emxPath)};
	(*dataRef)->open();

	ExpressionMatrix* matrix {dynamic_cast<ExpressionMatrix*>(&((*dataRef)->data()))};
	std::unique_ptr<float> expressions {matrix->dumpRawData()};

	// verify expression data
	float error = 0;

	for ( int i = 0; i < testExpressions.size(); ++i )
	{
		error += fabs(testExpressions[i] - expressions.get()[i]);
	}

	error /= testExpressions.size();

	Q_ASSERT(error < 1e-3);
}
