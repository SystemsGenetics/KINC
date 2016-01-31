/*
 * See GPL.txt for copyright information.
 *
 * Author: Joshua Burns
 *
 */
#ifndef CLDEVICE_H
#define CLDEVICE_H
#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#include <string>
#include "exception.h"



/// @brief Represents a single OpenCL device.
///
/// Represents a single OpenCL device that can be queried for information about
/// the device along with pointers to the OpenCL platform and device interface.
/// The "address" of these devices are represented by two numbers, the first
/// being the platform the device is part of and the second being which device
/// in the platform.
///
/// @author Josh Burns
/// @date 22 Jan 2016
class CLDevice
{
public:
   // *
   // * EXCEPTIONS
   // *
   struct OpenCLError;
   // *
   // * DECLERATIONS
   // *
   friend class CLDevList;
   // *
   // * ENUMERATIONS
   // *
   /// Defines all possible types of information that can be queried about an
   /// OpenCL device.
   enum CLInfo
   {
      ident, ///< The two identifying numbers of device.
      name, ///< The name of the device.
      type, ///< The OpenCL device type.
      online, ///< Is the device online or not.
      unified_mem, ///< Does the device have unified memory.
      addr_space, ///< The address space in bits.
      clock, ///< The maximum clock frequency.
      compute_units, ///< The number of compute units.
      work_size, ///< The maximum workgroup size.
      global_mem, ///< Size of global memory.
      local_mem, ///< Size of local memory.
      extensions, ///< Extensions this device supports.
      float_hp, ///< Floating point operations supported for half precision.
      float_sp, ///< Floating point operations supported for single precision.
      float_dp ///< Floating point operations supported for double precision.
   };
   // *
   // * FUNCTIONS
   // *
   cl::Platform& platform();
   cl::Device& device();
   std::string info(CLInfo) const;
   // *
   // * OPERATORS
   // *
   bool operator==(const CLDevice&);
private:
   // *
   // * BASIC METHODS
   // *
   CLDevice(int,int,const cl::Platform&,const cl::Device&);
   // *
   // * FUNCTIONS
   // *
   inline void yes_no(std::ostringstream&,bool) const;
   void size_it(std::ostringstream&,long) const;
   // *
   // * CONSTANTS
   // *
   static constexpr int _fPrecision {2};
   // *
   // * VARIABLES
   // *
   /// The increment into list of Platforms this device belongs to.
   int _pinc;
   /// The increment into list of devices of given platform that corresponds to
   /// this device.
   int _dinc;
   /// OpenCL platform of device.
   cl::Platform _platform;
   /// OpenCL object of actual device.
   cl::Device _device;
};



/// Exception is thrown when an OpenCL error occurs.
struct CLDevice::OpenCLError : public ::OpenCLError
{
   using ::OpenCLError::OpenCLError;
};



#endif
