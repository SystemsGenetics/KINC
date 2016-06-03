#ifndef CLPROGRAM_H
#define CLPROGRAM_H
#include <CL/cl.h>
#include <string>
#include <vector>
#include "exception.h"
#include "clkernel.h"



class CLProgram
{
public:
   ACE_EXCEPTION(CLProgram,NoSuchFile)
   ACE_EXCEPTION(CLProgram,NotCompiled)
   ACE_EXCEPTION(CLProgram,NotInitialized)
   OPENCL_EXCEPTION(CannotBind,clCreateProgramWithSource)
   OPENCL_EXCEPTION(BuildInfoFail,clGetProgramBuildInfo)
   OPENCL_EXCEPTION(CannotFindKern,clCreateKernel)
   CLProgram() = default;
   ~CLProgram();
   CLProgram(const CLProgram&) = delete;
   CLProgram& operator=(const CLProgram&) = delete;
   CLProgram(CLProgram&&) = delete;
   CLProgram& operator=(CLProgram&&) = delete;
protected:
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
