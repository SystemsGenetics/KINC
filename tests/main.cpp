#include "analyticfactory.h"
#include "datafactory.h"
#include "testclustermatrix.h"
#include "testcorrelationmatrix.h"
#include "testexpressionmatrix.h"



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

	ASSERT_TEST(new TestClusterMatrix);
	ASSERT_TEST(new TestCorrelationMatrix);
	ASSERT_TEST(new TestExpressionMatrix);

	return status;
}
