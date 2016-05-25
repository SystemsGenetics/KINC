#ifndef CLCONTEXT_H
#define CLCONTEXT_H
#include <CL/cl.h>
#include "cldevice.h"
#include "clprogram.h"
#include "clcommandqueue.h"
#include "clkernel.h"
#include "clbuffer.h"
#include "clevent.h"



class CLContext: public CLProgram, public CLCommandQueue
{
public:
   OPENCL_EXCEPTION(CreateContext,clCreateContext)
   ACE_EXCEPTION(CLContext,NotInitialized)
   CLContext() = default;
   ~CLContext();
   void init_cl(CLDevice& dev);
   bool is_initd();
protected:
   template<class T> CLBuffer<T> buffer(int);
private:
   bool _initd {false};
   cl_context _cid;
};



template<class T> CLBuffer<T> CLContext::buffer(int size)
{
   assert<NotInitialized>(_initd,__FILE__,__LINE__);
   return CLBuffer<T>(_cid,size);
}



#endif
