#ifndef ANALYTICPLUGIN_H
#define ANALYTICPLUGIN_H
#include "analytic.h"
#include "clprogram.h"
#include "clcommandqueue.h"
#include "clkernel.h"
#include "clbuffer.h"
#include "clevent.h"



class AnalyticPlugin : public Analytic, public CLProgram, CLCommandQueue
{
public:
   OPENCL_EXCEPTION(CreateContext,clCreateContext)
   ACE_EXCEPTION(AnalyticPlugin,NotInitialized)
   AnalyticPlugin() = default;
   ~AnalyticPlugin();
   void execute(GetOpts& ops, Terminal& tm);
   void init_cl(CLDevice& dev);
protected:
   template<class T> CLBuffer<T> buffer(int);
private:
   bool _isCL {false};
   cl_context _cid;
};



template<class T> CLBuffer<T> AnalyticPlugin::buffer(int size)
{
   assert<NotInitialized>(_isCL,__FILE__,__LINE__);
   return CLBuffer<T>(_cid,size);
}



#endif
