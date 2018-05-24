#include <ace/core/core.h>
#include <ace/core/ace_analytic_single.h>
#include <ace/core/ace_dataobject.h>

#include "testimportexpressionmatrix.h"
#include "analyticfactory.h"
#include "datafactory.h"
#include "importexpressionmatrix_input.h"



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

	// write expression data to text file
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

	// create analytic manager
	auto abstractManager = Ace::Analytic::AbstractManager::makeManager(AnalyticFactory::ImportExpressionMatrixType, 0, 1);
	auto manager = qobject_cast<Ace::Analytic::Single*>(abstractManager.release());
	manager->set(ImportExpressionMatrix::Input::InputFile, txtPath);
	manager->set(ImportExpressionMatrix::Input::OutputData, emxPath);
	manager->set(ImportExpressionMatrix::Input::NoSampleToken, noSampleToken);

	// run analytic
	manager->initialize();

	// TODO: wait for analytic to finish properly

	// read expression data from file
	std::unique_ptr<Ace::DataObject> dataRef {new Ace::DataObject(emxPath)};
	ExpressionMatrix* matrix {qobject_cast<ExpressionMatrix*>(dataRef->data())};
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
