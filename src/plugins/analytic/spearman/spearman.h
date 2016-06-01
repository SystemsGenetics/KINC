#ifndef SPEARMAN_H
#define SPEARMAN_H
#include "../../../analyticplugin.h"
#include "../../data/ematrix/ematrix.h"
#include "../../data/cmatrix/cmatrix.h"



class spearman : public AnalyticPlugin
{
public:
   ANALYTIC_EXCEPTION(spearman,InvalidInputType)
   ANALYTIC_EXCEPTION(spearman,TooManyOutputs)
   ANALYTIC_EXCEPTION(spearman,NotEnoughInputs)
   ANALYTIC_EXCEPTION(spearman,NoOutput)
   ANALYTIC_EXCEPTION(spearman,DifferentSampleSizes)
   void input(DataPlugin*) override final;
   void output(DataPlugin*) override final;
protected:
   void execute_cl(GetOpts&,Terminal&) override final;
   void execute_pn(GetOpts&,Terminal&) override final;
private:
   std::vector<ematrix*> _in;
   cmatrix* _out {nullptr};
};



#endif
