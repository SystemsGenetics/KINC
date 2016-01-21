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



class CLDevList;



/// @brief Represents a single OpenCL device.
///
/// Represents a single OpenCL device that can be queried for information about
/// the device along with pointers to the OpenCL platform and device interface.
/// The "address" of these devices are represented by two numbers, the first
/// being the platform the device is part of and the second being which device
/// in the platform.
class CLDevice
{
   friend class CLDevList;
public:
   // ****************************** Enumerations ***************************
   enum CLInfo {ident,name,type,online,unified_mem,addr_space,clock,
                compute_units,work_size,global_mem,local_mem,extensions,
                float_hp,float_sp,float_dp};
   // ****************************** Functions ******************************
   cl::Platform& platform();
   cl::Device& device();
   std::string info(CLInfo) const;
   // ****************************** Operators ******************************
   bool operator==(const CLDevice&);
private:
   // ****************************** Basic Methods **************************
   CLDevice(int,int,const cl::Platform&,const cl::Device&);
   // ****************************** Functions ******************************
   inline void yes_no(std::ostringstream&,bool) const;
   void size_it(std::ostringstream&,long) const;
   // ****************************** Constants ******************************
   static constexpr int _fPrecision {2};
   // ****************************** Variables ******************************
   int _pinc;
   int _dinc;
   cl::Platform _platform;
   cl::Device _device;
};



#endif
