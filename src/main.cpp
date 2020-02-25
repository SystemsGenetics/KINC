#if(GUI == 0)
#include <ace/cli/eapplication.h>
#else
#include <ace/gui/eapplication.h>
#include "gui/customizer.h"
#endif
#include "core/analyticfactory.h"
#include "core/datafactory.h"



using namespace std;



int main(int argc, char *argv[])
{
#if(GUI == 1)
    EAbstractCustomizer::setInstance(new Customizer);
#endif

    EApplication application(
        "SystemsGenetics"
        ,"kinc"
        ,KINC_MAJOR_VERSION
        ,KINC_MINOR_VERSION
        ,KINC_REVISION
        ,unique_ptr<DataFactory>(new DataFactory)
        ,unique_ptr<AnalyticFactory>(new AnalyticFactory)
        ,argc
        ,argv);
   return application.exec();
}
