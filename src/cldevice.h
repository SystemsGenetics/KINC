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
#include <vector>
#include <string>
#include <sstream>



enum class CLInfo { ident,name,type,online,unified_mem,addr_space,clock,
                    compute_units,work_size,global_mem,local_mem,extensions,
                    float_hp,float_sp,float_dp };



class CLDevice
/*
 * Represents a single OpenCL device that can be queried for information about
 * the device along with pointers to the OpenCL platform and device interface.
 * The "address" of these devices are represented by two numbers, the first
 * being the platform the device is part of and the second being which device in
 * the platform. Static functions and an iterator class is provided for
 * iterating through all available OpenCL devices.
 */
{
public:
   class Iterator;
   // ****************************** Functions ******************************
   static Iterator begin();
   static Iterator end();
   static bool exist(int,int);
   CLDevice() = delete;
   CLDevice(int,int);
   cl::Platform& platform();
   cl::Device& device();
   std::string info(CLInfo) const;
   // ****************************** Operators ******************************
   bool operator==(const CLDevice&);
private:
   // ****************************** Constants ******************************
   static constexpr int _fPrecision {2};
   // ****************************** Functions ******************************
   CLDevice(int,int,cl::Platform,cl::Device);
   inline void yesno(std::ostringstream&,bool) const;
   void sizeit(std::ostringstream&,long) const;
   // ****************************** Variables ******************************
   int _pinc;
   int _dinc;
   cl::Platform _platform;
   cl::Device _device;
};



class CLDevice::Iterator
/*
 * This iterates through all available OpenCL devices.
 */
{
   friend class CLDevice;
public:
   // ****************************** Functions ******************************
   static bool exist(int,int);
   Iterator() = delete;
   // ****************************** Operators ******************************
   CLDevice& operator*();
   void operator++();
   void operator++(int);
   bool operator!=(const Iterator&);
   bool operator==(const Iterator&);
private:
   // ****************************** Functions ******************************
   static void check_compiled();
   Iterator(bool);
   Iterator(int,int);
   // ****************************** Variables ******************************
   static bool _has_compiled;
   static std::vector<std::vector<CLDevice>> _devices;
   int _p;
   int _d;
};



#endif
