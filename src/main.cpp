#if(GUI == 0)
#include <ace/console/eapplication.h>
#else
#include <ace/gui/eapplication.h>
#endif
#include "core/analyticfactory.h"
#include "core/datafactory.h"



using namespace std;






int main(int argc, char *argv[])
{
   EApplication application("SystemsGenetics"
                            ,"kinc"
                            ,MAJOR_VERSION
                            ,MINOR_VERSION
                            ,REVISION
                            ,unique_ptr<DataFactory>(new DataFactory)
                            ,unique_ptr<AnalyticFactory>(new AnalyticFactory)
                            ,argc
                            ,argv);
   return application.exec();
}
