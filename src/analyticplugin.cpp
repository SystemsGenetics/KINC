#include "analyticplugin.h"



AnalyticPlugin::~AnalyticPlugin()
{
   if (_isCL)
   {
      clReleaseContext(_cid);
   }
}



void AnalyticPlugin::execute(GetOpts& ops, Terminal& tm)
{
   if (_isCL)
   {
      execute_cl(ops,tm);
   }
   else
   {
      execute_pn(ops,tm);
   }
}



void AnalyticPlugin::init_cl(CLDevice& dev)
{
   cl_context_properties props[] = {
      CL_CONTEXT_PLATFORM, (cl_context_properties)dev.platform(), 0 };
   cl_device_id device {dev.device()};
   cl_int err;
   _cid = clCreateContext(props,1,&device,NULL,NULL,&err);
   assert<CreateContext>(err==CL_SUCCESS,__FILE__,__LINE__,err);
   CLProgram::init(_cid,dev.device());
   _isCL = true;
}
