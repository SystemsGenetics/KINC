#include <ace/core/core.h>
#include <ace/core/ace_analytic_single.h>
#include <ace/core/ace_dataobject.h>

#include "testimportexpressionmatrix.h"
#include "../core/analyticfactory.h"
#include "../core/datafactory.h"
#include "../core/importexpressionmatrix_input.h"



void TestImportExpressionMatrix::test()
{
	// create random expression data
	int numGenes = 10;
	int numSamples = 5;
	std::vector<float> testExpressions(numGenes * numSamples);

	for ( size_t i = 0; i < testExpressions.size(); ++i )
	{
		testExpressions[i] = -10.0 + 20.0 * rand() / (1 << 31);
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

	// create raw text file
	QFile file(txtPath);
	QVERIFY(file.open(QIODevice::ReadWrite));

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
				stream << "\t" << nanToken;
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
	manager->set(ImportExpressionMatrix::Input::NANToken, nanToken);

	// run analytic
	manager->initialize();

	// TODO: wait for analytic to finish properly

	// read expression data from file
	std::unique_ptr<Ace::DataObject> dataRef {new Ace::DataObject(emxPath)};
	ExpressionMatrix* matrix {dataRef->data()->cast<ExpressionMatrix>()};
	std::vector<float> expressions {matrix->dumpRawData()};

	// verify expression data
	float error = 0;

	for ( size_t i = 0; i < testExpressions.size(); ++i )
	{
		error += fabs(testExpressions[i] - expressions[i]);
	}

	error /= testExpressions.size();

	QVERIFY(error < 1e-3);
}
