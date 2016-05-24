#ifndef ANALYTICPLUGIN_H
#define ANALYTICPLUGIN_H
#include "analytic.h"
#include "clprogram.h"
#include "clkernel.h"
#include "clbuffer.h"



class AnalyticPlugin : public Analytic, public CLProgram
{
public:
   OPENCL_EXCEPTION(CreateContext,clCreateContext)
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
   return CLBuffer<T>(_cid,size);
}



#endif
