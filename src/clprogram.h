#ifndef CLPROGRAM_H
#define CLPROGRAM_H
#include <CL/cl.h>
#include <string>
#include "cldevice.h"



class CLKernel;



class CLProgram
{
public:
   void init(CLDevice& dev) {}
   using string = std::string;
   bool add_source(const string& input, bool file = false);
   bool compile();
   string log();
   CLKernel mkernel(const string& name);
private:
   bool _initd {false};
};



#endif
