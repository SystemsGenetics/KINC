#ifdef ACE_GUI
   #include <ace/gui/Application.h>
#else
   #include <ace/console/Application.h>
#endif

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
