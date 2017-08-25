#include <ace/gui/Application.h>
#include "analyticfactory.h"
#include "datafactory.h"



using namespace std;



int main(int argc, char *argv[])
{
   unique_ptr<EAbstractAnalyticFactory> analyticFactory(new AnalyticFactory);
   unique_ptr<EAbstractDataFactory> dataFactory(new DataFactory);
   EAbstractAnalyticFactory::setInstance(move(analyticFactory));
   EAbstractDataFactory::setInstance(move(dataFactory));
   EApplication a(argc,argv,"KINC","kinc");
   return a.exec();
}
