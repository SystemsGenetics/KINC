#ifndef CLPROGRAM_H
#define CLPROGRAM_H
#include <CL/cl.h>
#include <string>
#include <vector>
#include "exception.h"



class CLKernel;



class CLProgram
{
public:
   ACE_EXCEPTION(CLProgram,NoSuchFile)
   OPENCL_EXCEPTION(CannotBind,clCreateProgramWithSource)
   OPENCL_EXCEPTION(BuildInfoFail,clGetProgramBuildInfo)
   CLProgram() = default;
   ~CLProgram();
   void init(cl_context,cl_device_id);
   using string = std::string;
   void add_source(const string& input, bool file = false);
   bool compile(const string&);
   string log();
   CLKernel mkernel(const string& name);
private:
   bool _initd {false};
   bool _binded {false};
   bool _compiled {false};
   cl_program _id;
   cl_device_id _did;
   cl_context _cid;
   std::vector<string> _src;
};



#endif
