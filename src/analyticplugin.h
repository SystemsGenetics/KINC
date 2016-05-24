#ifndef ANALYTICPLUGIN_H
#define ANALYTICPLUGIN_H
#include "analytic.h"
#include "clprogram.h"



class AnalyticPlugin : public Analytic, public CLProgram
{
public:
   void execute(GetOpts& ops, Terminal& tm);
   void init_cl(CLDevice& dev);
private:
   bool _isCL {false};
};



#endif
