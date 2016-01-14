/*
 * See GPL.txt for copyright information.
 *
 * Author: Joshua Burns
 *
 */
#include "cldevice.h"
#include "exception.h"
#include <iostream>



// Has a list of all OpenCL devices been built.
bool CLDevice::Iterator::_has_compiled {false};
// List of all OpenCL devices.
std::vector<std::vector<CLDevice>> CLDevice::Iterator::_devices {};



CLDevice::Iterator CLDevice::begin()
/*
 * Returns iterator at beginning of OpenCL device list.
 */
{
   return {false};
}



CLDevice::Iterator CLDevice::end()
/*
 * Returns iterator at one after end of OpenCL device list.
 */
{
   return {true};
}



bool CLDevice::exist(int p, int d)
/*
 * Verify if a certain OpenCL device exists.
 *
 * p: the increment into list of platforms.
 * d: the increment into list of devices for given platform.
 *
 * Returns true if OpenCL device does exist.
 */
{
   Iterator::exist(p,d);
}



CLDevice::CLDevice(int p, int d)
/*
 * Verifies that given OpenCL device exists, copies OpenCL platform and device
 * to new object.
 *
 * p: the increment into list of platforms.
 * d: the increment into list of devices for given platform.
 *
 * PRECONDITIONS:
 * 1. p must not exceed the last platform in list.
 * 2. d must not exceed the last device in given platform.
 */
   :_pinc(p),
   _dinc(d)
{
   Iterator i(p,d);
   _platform = (*i).platform();
   _device = (*i).device();
}



cl::Platform& CLDevice::platform()
/*
 * Return OpenCL platform of device.
 */
{
   return _platform;
}



cl::Device& CLDevice::device()
/*
 * Return OpenCL device of CLDevice object.
 */
{
   return _device;
}



std::string CLDevice::info(CLInfo which) const
/*
 * Returns string with requested information about OpenCL device.
 *
 * which: what type of data will be returned.
 */
{
   std::ostringstream buffer;
   try
   {
      switch (which)
      {
      case CLInfo::ident:
         buffer << _pinc << ":" << _dinc;
         break;
      case CLInfo::name:
         buffer << _device.getInfo<CL_DEVICE_NAME>();
         break;
      case CLInfo::type:
         switch (_device.getInfo<CL_DEVICE_TYPE>())
         {
         case CL_DEVICE_TYPE_CPU:
            buffer << "CPU";
            break;
         case CL_DEVICE_TYPE_GPU:
            buffer << "GPU";
            break;
         case CL_DEVICE_TYPE_ACCELERATOR:
            buffer << "ACCEL";
            break;
         case CL_DEVICE_TYPE_DEFAULT:
            buffer << "DEFAULT";
            break;
         default:
            buffer << "UNKNOWN";
            break;
         }
         break;
      case CLInfo::online:
         yesno(buffer,_device.getInfo<CL_DEVICE_AVAILABLE>()&&
                      _device.getInfo<CL_DEVICE_COMPILER_AVAILABLE>());
         break;
      case CLInfo::unified_mem:
         yesno(buffer,_device.getInfo<CL_DEVICE_HOST_UNIFIED_MEMORY>());
         break;
      case CLInfo::addr_space:
         buffer << cl_uint(_device.getInfo<CL_DEVICE_ADDRESS_BITS>());
         break;
      case CLInfo::clock:
         buffer << _device.getInfo<CL_DEVICE_MAX_CLOCK_FREQUENCY>();
         break;
      case CLInfo::compute_units:
         buffer << _device.getInfo<CL_DEVICE_MAX_COMPUTE_UNITS>();
         break;
      case CLInfo::work_size:
         buffer << _device.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
         break;
      case CLInfo::global_mem:
         sizeit(buffer,_device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>());
         break;
      case CLInfo::local_mem:
         sizeit(buffer,_device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>());
         break;
      case CLInfo::extensions:
         buffer << _device.getInfo<CL_DEVICE_EXTENSIONS>();
         break;
      }
   }
   catch (cl::Error e)
   {
      throw OpenCLError(__FILE__,__LINE__,e);
   }
   return buffer.str();
}



bool CLDevice::operator==(const CLDevice& cmp)
{
   return (_pinc==cmp._pinc&&_dinc==cmp._dinc);
}



CLDevice::CLDevice(int p, int d, cl::Platform pp, cl::Device dp):
   _pinc(p),
   _dinc(d),
   _platform(pp),
   _device(dp)
{}



inline void CLDevice::yesno(std::ostringstream& buffer, bool test) const
/*
 * Takes boolean value and outputs textual yes if boolean is true and no if
 * boolean is false.
 *
 * buffer: where textual yes/no will be written to.
 * test: boolean value that will be read.
 */
{
   const char* ans[] = {"yes","no"};
   buffer << ans[(test?0:1)];
}



void CLDevice::sizeit(std::ostringstream& buffer, long size) const
/*
 * Takes integer value as bytes and outputs formatted output to buffer with the
 * appropriate label. Maximum size is 1024 terrabytes.
 *
 * buffer: where formatted output is written.
 * size: size to be reported in bytes.
 */
{
   constexpr double kilobit = 1024.0f;
   constexpr int max = 4;
   const char* sizes[] = {"B","KB","MB","GB","TB"};
   int count = 0;
   double fsize = size;
   while (fsize>kilobit&&count<max)
   {
      fsize /= kilobit;
      count++;
   }
   buffer.setf(std::ios::fixed,std::ios::floatfield);
   buffer.precision(_fPrecision);
   buffer << fsize << sizes[count];
   buffer.unsetf(std::ios::floatfield);
}



bool CLDevice::Iterator::exist(int p, int d)
/*
 * Verify if a certain OpenCL device exists, making sure the OpenCL device list
 * has been generated beforehand.
 *
 * p: the increment into list of platforms.
 * d: the increment into list of devices for given platform.
 *
 * Returns true if OpenCL device does exist.
 */
{
   check_compiled();
   return (p<_devices.size()&&d<_devices[p].size());
}



CLDevice& CLDevice::Iterator::operator*()
{
   return _devices[_p][_d];
}



void CLDevice::Iterator::operator++()
{
   if (_devices.size()>0&&_d<_devices[_p].size())
   {
      _d++;
   }
   else
   {
      if (_p<_devices.size())
      {
         _d = 0;
         _p++;
      }
   }
}



void CLDevice::Iterator::operator++(int)
{
   ++(*this);
}



bool CLDevice::Iterator::operator!=(const Iterator& cmp)
{
   return (_p!=cmp._p||_d!=cmp._d);
}



bool CLDevice::Iterator::operator==(const Iterator& cmp)
{
   return (_p==cmp._p&&_d==cmp._d);
}



void CLDevice::Iterator::check_compiled()
/*
 * This static function will build the list of all available OpenCL devices.
 * The list is 2 dimensional, the first dimension the platform and the second
 * dimension the specific device of a platform. This list will only be built
 * once, the first time this static function is called.
 */
{
   if (!_has_compiled)
   {
      try
      {
         std::vector<cl::Platform> platforms;
         cl::Platform::get(&platforms);
         int i {0};
         for (auto p:platforms)
         {
            std::vector<cl::Device> devices;
            p.getDevices(CL_DEVICE_TYPE_ALL,&devices);
            _devices.push_back({});
            int j {0};
            for (auto d:devices)
            {
               _devices.back().push_back({i,j,p,d});
               j++;
            }
            i++;
         }
      }
      catch (cl::Error e)
      {
         throw OpenCLError(__FILE__,__LINE__,e);
      }

      _has_compiled = true;
   }
}



CLDevice::Iterator::Iterator(bool end)
/*
 * Checks to make sure OpenCL device list has been built first. Initializes
 * object to be at beginning of list or end depending on which one is specified.
 *
 * end: specifies what end of list to be initialized to, true if end and false
 *      if beginning.
 */
{
   check_compiled();
   if (end)
   {
      if (_devices.size()>0)
      {
         _p = _devices.size()-1;
         _d = _devices[_p].size();
      }
      else
      {
         _p = _d = 0;
      }
   }
   else
   {
      _p = _d = 0;
   }
}



CLDevice::Iterator::Iterator(int p, int d)
/*
 * Checks to make sure OpenCL device list has been built first. Initializes
 * object to point to specific device given.
 *
 * p: platform increment into list of devies.
 * d: device increment into list of given platform.
 *
 * PRECONDITIONS:
 * 1. p must not exceed the last platform in list.
 * 2. d must not exceed the last device in given platform.
 */
   :_p(p),
   _d(d)
{
   bool cond;
   check_compiled();
   cond = p<_devices.size()&&d<_devices[p].size();
   OutOfRange::assert(cond,__FILE__,__LINE__);
}
