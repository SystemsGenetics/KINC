#ifdef ACE_GUI
   #include <ace/gui/eapplication.h>
#else
   #include <ace/console/eapplication.h>
#endif

#include "analyticfactory.h"
#include "datafactory.h"



using namespace std;






int main(int argc, char *argv[])
{
   EApplication app(
      "", "kinc",
      3, 1, 2,
      unique_ptr<DataFactory>(new DataFactory),
      unique_ptr<AnalyticFactory>(new AnalyticFactory),
      argc, argv
   );

   return app.exec();
}
