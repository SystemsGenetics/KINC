#ifndef CLPROGRAM_H
#define CLPROGRAM_H
#include <CL/cl.h>
#include <string>



class CLKernel;



class CLProgram
{
public:
   void init(cl_context cid, cl_device_id did) {}
   using string = std::string;
   bool add_source(const string& input, bool file = false);
   bool compile();
   string log();
   CLKernel mkernel(const string& name);
private:
   bool _initd {false};
   bool _compiled {false};
   cl_program _id;
   cl_device_id _did;
};



#endif
