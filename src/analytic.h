#ifndef ANALYTIC_H
#define ANALYTIC_H
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#include <string>
#include "dataplugin.h"
#include "terminal.h"



class Analytic
{
public:
   virtual ~Analytic() = default;
   virtual void input(DataPlugin*) = 0;
   virtual void output(DataPlugin*) = 0;
   virtual void execute(GetOpts&,Terminal&,cl::Device*) = 0;
};



#endif
