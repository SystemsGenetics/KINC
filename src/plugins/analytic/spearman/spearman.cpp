#include "spearman.h"



void spearman::input(DataPlugin* input)
{
   bool cond {input->type()==std::string("ematrix")};
   assert<InvalidInputType>(cond,__FILE__,__LINE__);
   _in.push_back(input);
}



void spearman::output(DataPlugin* output)
{
   assert<TooManyOutputs>(_out==nullptr,__FILE__,__LINE__);
   _out = output;
}



void spearman::execute_cl(GetOpts& ops, Terminal& tm)
{
   ;
}



void spearman::execute_pn(GetOpts& ops, Terminal& tm)
{
   ;
}
