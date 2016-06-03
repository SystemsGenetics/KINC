/*
 * See GPL.txt for copyright information.
 *
 * Author: Joshua Burns
 *
 */
#ifndef CLDEVLIST_H
#define CLDEVLIST_H
#include <vector>
#include "exception.h"



class CLDevice;



/// @brief Container list of all OpenCL devices.
///
/// List of all available OpenCL devices. The list is constructed upon
/// instantiation of the object and can be recompiled by called the refresh()
/// function.
///
/// @author Josh Burns
/// @date 22 January 2016
class CLDevList
{
public:
   // *
   // * EXCEPTIONS
   // *
   OPENCL_EXCEPTION(PlatformErr,clGetPlatformIDs)
   OPENCL_EXCEPTION(DeviceErr,clGetDeviceIDs)
   // *
   // * DECLERATIONS
   // *
   class Iterator;
   // *
   // * BASIC METHODS
   // *
   CLDevList();
   // *
   // * FUNCTIONS
   // *
   Iterator begin();
   Iterator end();
   void refresh();
   bool exist(int,int);
   CLDevice& at(int,int);
private:
   // *
   // * FUNCTIONS
   // *
   void build();
   // *
   // * VARIABLES
   // *
   /// 2 dimensional list storing all OpenCL devices. The first vector
   /// corresponds to the list of all platforms and the second is the list of
   /// all devices for each platform.
   std::vector<std::vector<CLDevice>> _list;
};



/// Iterater for OpenCL device list.
///
/// @warning CLDevList object this iterate references cannot change during
/// lifetime of an iterator.
///
/// @author Josh Burns
/// @date 22 January 2016
class CLDevList::Iterator
{
public:
   // *
   // * DECLERATIONS
   // *
   friend class CLDevList;
   // *
   // * OPERATORS
   // *
   CLDevice& operator*();
   void operator++();
   void operator--();
   bool operator!=(const Iterator&);
private:
   // *
   // * BASIC METHODS
   // *
   Iterator(int,int,CLDevList*);
   // *
   // * VARIABLES
   // *
   /// Pointer to CLDevList object this iterator references.
   CLDevList* _devList;
   /// Current increment into list of platforms.
   int _pi;
   /// Current increment into list of devices of given platform.
   int _di;
};



#endif
