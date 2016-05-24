#ifndef CLKERNEL_H
#define CLKERNEL_H
#include <CL/cl.h>
#include "exception.h"



class CLKernel
{
   OPENCL_EXCEPTION(CannotSetArg,clSetKernelArg)
   OPENCL_EXCEPTION(CannotGetInfo,clGetKernelWorkGroupInfo)
   ACE_EXCEPTION(CLKernel,TooManyDims)
   ACE_EXCEPTION(CLKernel,DimOutOfRange)
   ACE_EXCEPTION(CLKernel,NotAlive)
   constexpr static int _maxDims {16};
   friend class CLCommadQueue;
   friend class CLProgram;
   ~CLKernel();
   CLKernel(const CLKernel&);
   CLKernel& operator=(const CLKernel&);
   CLKernel(CLKernel&&);
   CLKernel& operator=(CLKernel&&);
   template<class T> void set_arg(cl_uint,T);
   template<class... Args> void set_args(Args...);
   void set_swarm_dims(cl_uint);
   void set_swarm_size(int,cl_uint,cl_uint);
   size_t get_wg_size();
   size_t get_wg_multiple();
private:
   CLKernel(cl_kernel,cl_device_id);
   template<class T> void set_args_int(T,int);
   template<class T, class... Args> void set_args_int(T,int,Args...);
   bool _isAlive {false};
   cl_device_id _did;
   cl_kernel _id;
   cl_uint _dims {1};
   size_t _gSizes[_maxDims] {1};
   size_t _lSizes[_maxDims] {1};
};



template<class T> void CLKernel::set_arg(cl_uint index, T arg)
{
   assert<NotAlive>(_isAlive,__FILE__,__LINE__);
   cl_int err = clSetKernelArg(_id,index,sizeof(T),&arg);
   assert<CannotSetArg>(err==CL_SUCCESS,__FILE__,__LINE__,err);
}



template<class... Args> void CLKernel::set_args(Args... args)
{
   assert<NotAlive>(_isAlive,__FILE__,__LINE__);
   set_args_int(args...,0);
}



template<class T> void CLKernel::set_args_int(T arg, int ind)
{
   set_arg(ind,arg);
}



template<class T, class... Args>
void CLKernel::set_args_int(T arg, int ind, Args... args)
{
   set_arg(ind,arg);
   set_args(args...,ind+1);
}



#endif
