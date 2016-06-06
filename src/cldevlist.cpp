#include <CL/cl.h>
#include <memory>
#include "cldevlist.h"
#include "cldevice.h"



/// Compile list of available OpenCL devices.
CLDevList::CLDevList()
{
   build();
}



/// @brief Builds list of OpenCL devices for object.
///
/// This internal function constructs all available OpenCL devices by querying
/// OpenCL and constructing a 2 dimensional vector storing all devices.
///
/// @pre The internal vector list must be empty.
void CLDevList::build()
{
   cl_platform_id* platforms {nullptr};
   cl_device_id* devices {nullptr};
   try
   {
      cl_uint ptotal;
      cl_int err = clGetPlatformIDs(0,NULL,&ptotal);
      assert<PlatformErr>(err==CL_SUCCESS,__FILE__,__LINE__,err);
      platforms = new cl_platform_id[ptotal];
      err = clGetPlatformIDs(ptotal,platforms,NULL);
      assert<PlatformErr>(err==CL_SUCCESS,__FILE__,__LINE__,err);
      for (int i=0;i<ptotal;++i)
      {
         cl_uint dtotal;
         err = clGetDeviceIDs(platforms[i],CL_DEVICE_TYPE_ALL,0,NULL,&dtotal);
         assert<DeviceErr>(err==CL_SUCCESS,__FILE__,__LINE__,err);
         devices = new cl_device_id[dtotal];
         err = clGetDeviceIDs(platforms[i],CL_DEVICE_TYPE_ALL,dtotal,devices,
                              NULL);
         assert<DeviceErr>(err==CL_SUCCESS,__FILE__,__LINE__,err);
         _list.push_back({});
         for (int j=0;j<dtotal;++j)
         {
            _list.back().emplace_back(i,j,platforms[i],devices[j]);
         }
         delete[] devices;
         devices = nullptr;
      }
      delete[] platforms;
   }
   catch (...)
   {
      _list.clear();
      if (platforms) delete[] platforms;
      if (devices) delete[] devices;
      throw;
   }
}



/// Return iterator at beginning of list.
CLDevList::Iterator CLDevList::begin()
{
   return {0,0,this};
}



/// Return iterator at one past end of list.
CLDevList::Iterator CLDevList::end()
{
   int pi = 0;
   int di = 0;
   if (_list.size()>0)
   {
      pi = static_cast<int>(_list.size()-1);
      di = static_cast<int>(_list.back().size());
   }
   return {pi,di,this};
}



/// @brief Refresh OpenCL device list.
///
/// Clears vector list of previously built devices and builds new list of OpenCL
/// devices.
void CLDevList::refresh()
{
   _list.clear();
   build();
}



/// Returns if specified OpenCL device exists.
///
/// @param p The increment into list of platforms.
/// @param d The increment into list of devices of given platform.
/// @return True if given OpenCL device exists.
bool CLDevList::exist(int p, int d)
{
   return (p<_list.size()&&d<_list[p].size());
}



/// Returns specified OpenCL device.
///
/// @param p The increment into list of platforms.
/// @param d The increment into list of devices of given platform.
/// @return OpenCL device.
///
/// @pre p must not exceed the last platform.
/// @pre d must not exceed the last device in given platform.
CLDevice& CLDevList::at(int p, int d)
{
   return _list.at(p).at(d);
}



/// Returns reference to current OpenCL device.
CLDevice& CLDevList::Iterator::operator*()
{
   return _devList->_list[_pi][_di];
}



/// @brief Increment to next OpenCL device in list.
///
/// Will move to next device in list of OpenCL devices unless iterator is
/// already at end of list in which case it stays the same.
void CLDevList::Iterator::operator++()
{
   if (_devList->_list.size()>0)
   {
      if (++_di>=_devList->_list[_pi].size())
      {
         _di = 0;
         if (++_pi>=_devList->_list.size())
         {
            _pi--;
            _di = _devList->_list.back().size();
         }
      }
   }
}



/// @brief Decrement to previous OpenCL device in list.
///
/// Will move to previous device in list of OpenCL devices unless iterator is
/// already at beginning of list in which case it stays the same.
void CLDevList::Iterator::operator--()
{
   if (--_di==-1)
   {
      _di = 0;
      if (--_pi==-1)
      {
         _pi = 0;
      }
   }
}



/// Defines inequality between two iterators.
bool CLDevList::Iterator::operator!=(const Iterator& cmp)
{
   return (_pi!=cmp._pi||_di!=cmp._di);
}



/// Initializes iterator with given parameters.
///
/// @param p The increment into list of platforms.
/// @param d The increment into list of devices of given platform.
/// @param devList Pointer to CLDevList object which instantiated iterator.
///
/// @pre p must not exceed the last platform.
/// @pre d must not exceed the last device in given platform.
CLDevList::Iterator::Iterator(int p, int d, CLDevList* devList):
   _pi(p),
   _di(d),
   _devList(devList)
{}
