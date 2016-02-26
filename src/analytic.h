#ifndef ANALYTIC_H
#define ANALYTIC_H
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#include <string>
#include <memory>
#include "dataplugin.h"
#include "terminal.h"



class Analytic
{
public:
   using wptr = std::weak_ptr<DataPlugin>;
   virtual ~Analytic() = default;
   virtual void input(wptr) = 0;
   virtual void output(wptr) = 0;
   virtual void execute(GetOpts&,Terminal&,cl::Device*) = 0;
};



#endif
