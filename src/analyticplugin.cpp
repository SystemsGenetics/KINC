#include "analyticplugin.h"



void AnalyticPlugin::execute(GetOpts& ops, Terminal& tm)
{
   if (CLContext::is_initd())
   {
      execute_cl(ops,tm);
   }
   else
   {
      execute_pn(ops,tm);
   }
}
