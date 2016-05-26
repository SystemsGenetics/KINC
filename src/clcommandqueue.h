#ifndef CLCOMMANDQUEUE_H
#define CLCOMMANDQUEUE_H
#include <CL/cl.h>
#include "clkernel.h"
#include "clbuffer.h"
#include "clevent.h"



class CLCommandQueue
{
public:
   OPENCL_EXCEPTION(CannotCreate,clCreateCommandQueue)
   OPENCL_EXCEPTION(CannotEnqueueRB,clEnqueueReadBuffer)
   OPENCL_EXCEPTION(CannotEnqueueWB,clEnqueueWriteBuffer)
   OPENCL_EXCEPTION(CannotAddTask,clEnqueueTask)
   OPENCL_EXCEPTION(CannotAddSwarm,clEnqueueNDRangeKernel)
   ACE_EXCEPTION(CLCommandQueue,NotInitialized)
   ACE_EXCEPTION(CLCommandQueue,DeadBufferUsed)
   ACE_EXCEPTION(CLCommandQueue,DeadKernelUsed)
   CLCommandQueue() = default;
   ~CLCommandQueue();
   CLCommandQueue(const CLCommandQueue&) = delete;
   CLCommandQueue& operator=(const CLCommandQueue&) = delete;
   CLCommandQueue(CLCommandQueue&&) = delete;
   CLCommandQueue& operator=(CLCommandQueue&&) = delete;
protected:
   void init(cl_context,cl_device_id);
   template<class T> CLEvent read_buffer(CLBuffer<T>&);
   template<class T> CLEvent write_buffer(CLBuffer<T>&);
   CLEvent add_task(CLKernel&);
   CLEvent add_swarm(CLKernel&);
private:
   bool _initd {false};
   cl_command_queue _id;
};



template<class T> CLEvent CLCommandQueue::read_buffer(CLBuffer<T>& buffer)
{
   assert<NotInitialized>(_initd,__FILE__,__LINE__);
   assert<DeadBufferUsed>(buffer._hostPtr,__FILE__,__LINE__);
   cl_event ret;
   cl_int err = clEnqueueReadBuffer(_id,buffer._id,CL_FALSE,0,
                                    buffer._size*sizeof(T),buffer._hostPtr,0,
                                    NULL,&ret);
   assert<CannotEnqueueRB>(err==CL_SUCCESS,__FILE__,__LINE__,err);
   return CLEvent(ret);
}



template<class T> CLEvent CLCommandQueue::write_buffer(CLBuffer<T>& buffer)
{
   assert<NotInitialized>(_initd,__FILE__,__LINE__);
   assert<DeadBufferUsed>(buffer._hostPtr,__FILE__,__LINE__);
   cl_event ret;
   cl_int err = clEnqueueWriteBuffer(_id,buffer._id,CL_FALSE,0,
                                     buffer._size*sizeof(T),buffer._hostPtr,0,
                                     NULL,&ret);
   assert<CannotEnqueueWB>(err==CL_SUCCESS,__FILE__,__LINE__,err);
   return CLEvent(ret);
}



#endif
