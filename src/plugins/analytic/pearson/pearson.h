#ifndef PEARSON_H
#define PEARSON_H
#include "../../../analyticplugin.h"



class pearson : public AnalyticPlugin
{
public:
   ANALYTIC_EXCEPTION(pearson,InvalidInputType)
   ANALYTIC_EXCEPTION(pearson,TooManyOutputs)
   void input(DataPlugin*) override final;
   void output(DataPlugin*) override final;
protected:
   void execute_cl(GetOpts&,Terminal&) override final;
   void execute_pn(GetOpts&,Terminal&) override final;
private:
   std::vector<DataPlugin*> _in;
   DataPlugin* _out {nullptr};
};



#endif
