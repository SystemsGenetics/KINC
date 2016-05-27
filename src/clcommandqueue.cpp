#include "clcommandqueue.h"



CLCommandQueue::~CLCommandQueue()
{
   if (_initd)
   {
      clReleaseCommandQueue(_id);
   }
}



void CLCommandQueue::init(cl_context cid, cl_device_id did)
{
   cl_int err;
   _id = clCreateCommandQueue(cid,did,0x0,&err);
   assert<CannotCreate>(err==CL_SUCCESS,__FILE__,__LINE__,err);
   _initd = true;
}



CLEvent CLCommandQueue::add_task(CLKernel& kernel)
{
   assert<NotInitialized>(_initd,__FILE__,__LINE__);
   assert<DeadKernelUsed>(kernel._isAlive,__FILE__,__LINE__);
   cl_event ret;
   cl_int err;
   err = clEnqueueTask(_id,kernel._id,0,NULL,&ret);
   assert<CannotAddTask>(err==CL_SUCCESS,__FILE__,__LINE__,err);
   return CLEvent(ret);
}



CLEvent CLCommandQueue::add_swarm(CLKernel& kernel)
{
   assert<NotInitialized>(_initd,__FILE__,__LINE__);
   assert<DeadKernelUsed>(kernel._isAlive,__FILE__,__LINE__);
   cl_event ret;
   const size_t offsets[kernel._dims] = {0};
   cl_int err;
   err = clEnqueueNDRangeKernel(_id,kernel._id,kernel._dims,offsets,
                                kernel._gSizes,kernel._lSizes,
                                0,NULL,&ret);
   assert<CannotAddSwarm>(err==CL_SUCCESS,__FILE__,__LINE__,err);
   return CLEvent(ret);
}
