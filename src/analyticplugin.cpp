#include "analyticplugin.h"



void AnalyticPlugin::execute(GetOpts& ops, Terminal& tm)
{
   if (_isCL)
   {
      execute_cl(ops,tm);
   }
   else
   {
      execute_pn(ops,tm);
   }
}



void AnalyticPlugin::init_cl(CLDevice& dev)
{
   CLProgram::init(dev);
   _isCL = true;
}
