#ifndef SPEARMAN_H
#define SPEARMAN_H
#include "../../../analytic.h"



class spearman : public Analytic
{
public:
   ALTC_EXCEPTION(pearson,InvalidInputType)
   ALTC_EXCEPTION(pearson,TooManyOutputs)
   void input(DataPlugin*) override final;
   void output(DataPlugin*) override final;
   void execute(GetOpts&,Terminal&,cl::Device*) override final;
   void execute(GetOpts&,Terminal&) override final;
private:
   std::vector<DataPlugin*> _in;
   DataPlugin* _out {nullptr};
};



#endif
