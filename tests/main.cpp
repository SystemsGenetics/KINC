#include "analyticfactory.h"
#include "datafactory.h"
#include "testclustermatrix.h"
#include "testcorrelationmatrix.h"
#include "testexportexpressionmatrix.h"
#include "testexpressionmatrix.h"
#include "testimportexpressionmatrix.h"



int main(int argc, char **argv)
{
	std::unique_ptr<EAbstractAnalyticFactory> analyticFactory(new AnalyticFactory);
   std::unique_ptr<EAbstractDataFactory> dataFactory(new DataFactory);
   EAbstractAnalyticFactory::setInstance(move(analyticFactory));
   EAbstractDataFactory::setInstance(move(dataFactory));

	int status {0};
	auto ASSERT_TEST = [&status, argc, argv](QObject* object)
	{
		status |= QTest::qExec(object, argc, argv);
		delete object;
	};

	try
	{
		ASSERT_TEST(new TestClusterMatrix);
		ASSERT_TEST(new TestCorrelationMatrix);
		ASSERT_TEST(new TestExportExpressionMatrix);
		ASSERT_TEST(new TestExpressionMatrix);
		ASSERT_TEST(new TestImportExpressionMatrix);
	}
	catch ( EException& e )
	{
		QTextStream stream(stdout);
		stream << QObject::tr("CRITICAL ERROR\n\n");
		stream << e.getTitle() << QObject::tr("\n\n");
		stream << e.getDetails() << QObject::tr("\n\n");
		stream << QObject::tr("File: %1\nLine: %2\nFunction: %3\n")
			.arg(e.getFile())
			.arg(e.getLine())
			.arg(e.getFunction());
	}

	return status;
}
