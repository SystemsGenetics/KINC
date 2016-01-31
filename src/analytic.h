#ifndef ANALYTIC_H
#define ANALYTIC_H
#include <string>



class Analytic
{
   virtual ~Analytic() = default;
   virtual bool input(const std::string&) = 0;
   virtual bool output(const std::string&) = 0;
   virtual void option(const std::string&,const std::string&) = 0;
   virtual bool execute() = 0;
};



#endif
