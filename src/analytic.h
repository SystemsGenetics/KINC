#ifndef ANALYTIC_H
#define ANALYTIC_H
#include <string>
#include "dataplugin.h"



class Analytic
{
   virtual ~Analytic() = default;
   virtual void input(DataPlugin*) = 0;
   virtual void output(DataPlugin*) = 0;
   virtual void option(const std::string&,const std::string&) = 0;
   virtual void execute() = 0;
};



#endif
