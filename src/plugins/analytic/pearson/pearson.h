#ifndef EMATRIX_H
#define EMATRIX_H
#include <fstream>
#include "../../../analytic.h"



class pearson : public Analytic
{
public:
   void input(DataPlugin*) override final {}
   void output(DataPlugin*) override final {}
   void execute(GetOpts&,Terminal&,cl::Device*) override final {}
};



#endif
