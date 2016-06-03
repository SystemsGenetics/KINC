#include "clcontext.h"



CLContext::~CLContext()
{
   if (_initd)
   {
      clReleaseContext(_cid);
   }
}



void CLContext::init_cl(CLDevice& dev)
{
   cl_context_properties props[] = {
      CL_CONTEXT_PLATFORM, (cl_context_properties)dev.platform(), 0 };
   cl_device_id device {dev.device()};
   cl_int err;
   _cid = clCreateContext(props,1,&device,NULL,NULL,&err);
   assert<CreateContext>(err==CL_SUCCESS,__FILE__,__LINE__,err);
   try
   {
      CLProgram::init(_cid,dev.device());
      CLCommandQueue::init(_cid,dev.device());
   }
   catch (...)
   {
      clReleaseContext(_cid);
      throw;
   }
   _initd = true;
}



bool CLContext::is_initd()
{
   return _initd;
}
