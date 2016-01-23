/*
 * See GPL.txt for copyright information.
 *
 * Author: Joshua Burns
 *
 */
#include <sstream>
#include "cldevice.h"
#include "exception.h"



/// Return OpenCL platform of device.
cl::Platform& CLDevice::platform()
{
   return _platform;
}



/// Return OpenCL device of CLDevice object.
cl::Device& CLDevice::device()
{
   return _device;
}



/// @brief Returns information about OpenCL device.
///
/// Returns string with requested information about OpenCL device. The different
/// types of information to be queired is defined as a class enumeration.
///
/// @param which Enumerated type specifying what type of information to query.
std::string CLDevice::info(CLInfo which) const
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
         yes_no(buffer,_device.getInfo<CL_DEVICE_AVAILABLE>()&&
                      _device.getInfo<CL_DEVICE_COMPILER_AVAILABLE>());
         break;
      case CLInfo::unified_mem:
         yes_no(buffer,_device.getInfo<CL_DEVICE_HOST_UNIFIED_MEMORY>());
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
         size_it(buffer,_device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>());
         break;
      case CLInfo::local_mem:
         size_it(buffer,_device.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>());
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



/// @brief Constructs new CLDevice based off info given.
///
/// This is the only constructor for this class which is private and should only
/// be used by the CLDevList class. The parameters given for the device are NOT
/// checked for validity.
///
/// @param p The increment into list of platforms.
/// @param d The increment into list of devices for given platform.
/// @param pl The OpenCL platform object for the device.
/// @param de The OpenCL device object for the device.
///
/// @pre p must not exceed the last platform.
/// @pre d must not exceed the last device in given platform.
CLDevice::CLDevice(int p, int d, const cl::Platform& pl, const cl::Device& de)
   :_pinc(p),
   _dinc(d),
   _platform(pl),
   _device(de)
{}



/// @brief Converts boolean value into string representation.
///
/// Takes boolean value and outputs textual yes if boolean is true and no if
/// boolean is false.
///
/// @param buffer Where textual yes/no will be written to.
/// @param test Value that will be read.
inline void CLDevice::yes_no(std::ostringstream& buffer, bool test) const
{
   const char* ans[] = {"yes","no"};
   buffer << ans[(test?0:1)];
}



void CLDevice::size_it(std::ostringstream& buffer, long size) const
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
