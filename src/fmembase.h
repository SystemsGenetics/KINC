#ifndef FILEMEM_H
#define FILEMEM_H
#include <string>
#include <cstdint>
#include <type_traits>
#include <utility>
#include <cstring>
#include "exception.h"



class LinuxFile;



class FileMem
{
public:
   // *
   // * EXCEPTIONS
   // *
   struct Exception;
   struct SystemError;
   struct InvalidFile;
   struct FileSegFault;
   struct OutOfMemory;
   // *
   // * DECLERATIONS
   // *
   enum class Sync { read,write };
   using Base = LinuxFile;
   using VPtr = uint64_t;
   template<int S> struct Object;
   struct DObject;
   template<class M,class T> class RawMap;
   template<class T> using Map = RawMap<LinuxFile,T>;
   template<class M> struct RawPtr
   {
      VPtr inc;
      M* fmem;
   };
   using Ptr = RawPtr<LinuxFile>;
   // *
   // * CONSTANTS
   // *
   /// Defines the value of a null pointer for file memory.
   constexpr static VPtr nullPtr {0xffffffffffffffffull};
   /// Identifying string at the beginning of any memory object file.
   constexpr const static char* _identString = "\0\15\41\102\104\101\124\0\0";
   constexpr static int _idLen = 9;
};



template<int S> struct FileMem::Object
{
   template<class T> inline T& get(int n);
   char bytes[S];
   constexpr static uint64_t size = S;
};



struct FileMem::DObject
{
   inline DObject(const DObject&);
   inline DObject(DObject&&);
   inline DObject& operator=(const DObject&);
   inline DObject& operator=(DObject&&);
   inline DObject(uint64_t);
   template<class T> inline T& get(int);
   inline ~DObject();
   char* bytes;
   uint64_t size;
};



template<class M,class T> class FileMem::RawMap
{
public:
   inline void operator=(VPtr);
   inline RawMap(RawPtr<M>);
   inline VPtr addr() const;
   inline void sync(Sync,uint64_t = 0);
   inline T& operator*();
   inline T* operator->();
private:
   T _val;
   RawPtr<M> _ptr;
};



template<int S> template<class T> inline T& FileMem::Object<S>::get(int n)
{
   return *reinterpret_cast<T*>(&bytes[n]);
}



inline FileMem::DObject::DObject(const DObject& obj):
   size(obj.size)
{
   bytes = new char[size];
   memcpy(bytes,obj.bytes,size);
}



inline FileMem::DObject::DObject(DObject&& tmp):
   size(tmp.size),
   bytes(tmp.bytes)
{
   tmp.bytes = nullptr;
   tmp.size = 0;
}



inline FileMem::DObject& FileMem::DObject::operator=(const DObject& obj)
{
   size = obj.size;
   if (bytes)
   {
      delete[] bytes;
   }
   bytes = new char[size];
   memcpy(bytes,obj.bytes,size);
}



inline FileMem::DObject& FileMem::DObject::operator=(DObject&& tmp)
{
   size = tmp.size;
   bytes = tmp.bytes;
   tmp.bytes = nullptr;
   tmp.size = 0;
}



inline FileMem::DObject::DObject(uint64_t dSize):
   size(dSize)
{
   bytes = new char[size];
}



inline FileMem::DObject::~DObject()
{
   if (bytes)
   {
      delete[] bytes;
   }
}



template<class T> inline T& FileMem::DObject::get(int n)
{
   return *reinterpret_cast<T*>(&bytes[n]);
}



template<class M,class T> inline void FileMem::RawMap<M,T>::operator=(VPtr vptr)
{
   _ptr.inc = vptr;
}



template<class M,class T> inline FileMem::RawMap<M,T>::RawMap(RawPtr<M> loc):
   _ptr(loc)
{}



template<class M,class T>
inline FileMem::VPtr FileMem::RawMap<M,T>::addr() const
{
   return _ptr.inc;
}



template<class M,class T> inline void FileMem::RawMap<M,T>::sync(Sync which,
                                                                 uint64_t inc)
{
   uint64_t offset = _ptr.inc + inc*_val.size;
   switch (which)
   {
   case Sync::read:
      _ptr.fmem->read(_val.bytes,offset,_val.size);
      break;
   case Sync::write:
      _ptr.fmem->write(_val.bytes,offset,_val.size);
      break;
   }
}



template<class M,class T> inline T& FileMem::RawMap<M,T>::operator*()
{
   return _val;
}



template<class M,class T> inline T* FileMem::RawMap<M,T>::operator->()
{
   return &_val;
}



/// Generic base exception class for all exceptions thrown in FileMem class.
struct FileMem::Exception : public ::Exception
{
   using ::Exception::Exception;
};

/// Exception that is thrown when a system error occurs.
struct FileMem::SystemError : public ::SystemError
{
   using ::SystemError::SystemError;
};

/// Exception that is thrown when a file that is not a valid file memory object
/// is opened.
struct FileMem::InvalidFile : public FileMem::Exception
{
   InvalidFile(const char* file, int line):
      Exception(file,line,"FileMem::InvalidFile")
   {}
};

/// Exception that is thrown when a read or write operation is attempted that
/// goes beyond the limits of the total file object's size.
struct FileMem::FileSegFault : public FileMem::Exception
{
   FileSegFault(const char* file, int line):
      Exception(file,line,"FileMem::FileSegFault")
   {}
};

/// Exception that is thrown when an allocation of new file memory is attempted
/// that is greater than the available amount of bytes that can be allocated.
struct FileMem::OutOfMemory : public FileMem::Exception
{
   OutOfMemory(const char* file, int line):
      Exception(file,line,"FileMem::OutOfMemory")
   {}
};



#endif
