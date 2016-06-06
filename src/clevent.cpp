#include "clevent.h"



CLEvent::~CLEvent()
{
   if (_hasEvent)
   {
      clReleaseEvent(_id);
   }
}



CLEvent::CLEvent(CLEvent&& move):
   _hasEvent(move._hasEvent),
   _id(move._id)
{
   move._hasEvent = false;
}



CLEvent& CLEvent::operator=(CLEvent&& move)
{
   if (_hasEvent)
   {
      clReleaseEvent(_id);
   }
   _hasEvent = move._hasEvent;
   _id = move._id;
   move._hasEvent = false;
}



void CLEvent::wait()
{
   if (_hasEvent)
   {
      cl_int err;
      err = clWaitForEvents(1,&_id);
      assert<CannotWait>(err==CL_SUCCESS,__FILE__,__LINE__,err);
   }
}



bool CLEvent::is_done()
{
   bool ret {true};
   if (_hasEvent)
   {
      cl_int status;
      cl_int err;
      err = clGetEventInfo(_id,CL_EVENT_COMMAND_EXECUTION_STATUS,sizeof(cl_int),
                           &status,NULL);
      assert<CannotGetInfo>(err==CL_SUCCESS,__FILE__,__LINE__,err);
      if (status<0)
      {
         throw ExecutionFail(__FILE__,__LINE__);
      }
      if (status!=CL_COMPLETE)
      {
         ret = false;
      }
   }
   return ret;
}



CLEvent::CLEvent(cl_event id):
   _hasEvent(true),
   _id(id)
{}
