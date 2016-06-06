#include "pearson.h"



void pearson::input(DataPlugin* input)
{
   bool cond {input->type()==std::string("ematrix")};
   assert<InvalidInputType>(cond,__FILE__,__LINE__);
   _in.push_back(input);
}



void pearson::output(DataPlugin* output)
{
   assert<TooManyOutputs>(_out==nullptr,__FILE__,__LINE__);
   _out = output;
}



void pearson::execute_cl(GetOpts& ops, Terminal& tm)
{
   ;
}



void pearson::execute_pn(GetOpts& ops, Terminal& tm)
{
   ;
}
