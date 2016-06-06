#ifndef CLEVENT_H
#define CLEVENT_H
#include <CL/cl.h>
#include "exception.h"



class CLEvent
{
public:
   OPENCL_EXCEPTION(CannotWait,clWaitForEvents)
   OPENCL_EXCEPTION(CannotGetInfo,clGetEventInfo)
   ACE_EXCEPTION(CLEvent,ExecutionFail)
   friend class CLCommandQueue;
   CLEvent() = default;
   ~CLEvent();
   CLEvent(const CLEvent&) = delete;
   CLEvent& operator=(const CLEvent&) = delete;
   CLEvent(CLEvent&&);
   CLEvent& operator=(CLEvent&&);
   void wait();
   bool is_done();
private:
   CLEvent(cl_event);
   bool _hasEvent {false};
   cl_event _id;
};



#endif
