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
   using Base = LinuxFile;
   using VPtr = uint64_t;
   template<int S> struct Object;
   struct DObject;
   template<class M,class T> class RawPtr;
   template<class T> using Ptr = RawPtr<LinuxFile,T>;
   // *
   // * CONSTANTS
   // *
   /// Defines the value of a null pointer for file memory.
   constexpr static VPtr nullPtr = 0xffffffffffffffffll;
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



template<class M,class T> class FileMem::RawPtr
{
public:
   static_assert(std::is_class<T>::value,
                 "FileMem::Ptr<T>, T is of an unsupported type.");
   RawPtr(const Ptr<T>&) = delete;//SADFACE... AND THIS :(((
   inline RawPtr(RawPtr<M,T>&&);
   RawPtr<M,T>& operator=(const Ptr<T>&) = delete; // SAD FACE, we will need this
   inline void operator=(RawPtr<M,T>&&);
   inline void operator=(VPtr);
   inline RawPtr(M&,VPtr = nullPtr,uint64_t = 1);
   inline VPtr addr() const;
   inline void raw(uint64_t);
   inline void save();
   inline T& operator*();
   inline T* operator->();
   inline T& operator[](uint64_t);
   void operator++();
   void operator--();
private:
   T _val;
   M* _base;
   VPtr _ptr;
   uint64_t _size;
   uint64_t _inc;
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



template<class M,class T>
inline FileMem::RawPtr<M,T>::RawPtr(FileMem::RawPtr<M,T>&& tmp):
   _val(std::move(tmp._val)),
   _base(tmp._base),
   _ptr(tmp._ptr),
   _size(tmp._size),
   _inc(tmp._inc)
{
   tmp._ptr = nullPtr;
}



template<class M,class T>
inline void FileMem::RawPtr<M,T>::operator=(RawPtr<M,T>&& tmp)
{
   _val = std::move(tmp._val);
   _base = tmp._base;
   _ptr = tmp._ptr;
   _size = tmp._size;
   _inc = tmp._inc;
   tmp._ptr = nullPtr;
}



template<class M,class T>
inline void FileMem::RawPtr<M,T>::operator=(VPtr ptr)
{
   _ptr = ptr;
   _inc = 0;
   _base->read(_val.bytes,_ptr,_val.size);
}



template<class M,class T> inline FileMem::RawPtr<M,T>::RawPtr(M& mem, VPtr loc,
                                                              uint64_t size):
   _base(&mem),
   _ptr(loc),
   _size(size),
   _inc(0)
{
   assert<FileSegFault>(size!=0,__FILE__,__LINE__);
   if (loc!=nullPtr)
   {
      _base->read(_val.bytes,_ptr,_val.size);
   }
   else
   {
      _ptr = _base->allocate(_val.size*_size);
   }
}



template<class M,class T>
inline FileMem::VPtr FileMem::RawPtr<M,T>::addr() const
{
   return _ptr;
}



template<class M,class T> inline void FileMem::RawPtr<M,T>::raw(uint64_t inc)
{
   assert<FileSegFault>(inc<_size,__FILE__,__LINE__);
   _inc = inc;
}



template<class M,class T> inline void FileMem::RawPtr<M,T>::save()
{
   _base->write(_val.bytes,_ptr + _inc*_val.size,_val.size);
}



template<class M,class T> inline T& FileMem::RawPtr<M,T>::operator*()
{
   return _val;
}



template<class M,class T> inline T* FileMem::RawPtr<M,T>::operator->()
{
   return &_val;
}



template<class M,class T>
inline T& FileMem::RawPtr<M,T>::operator[](uint64_t inc)
{
   assert<FileSegFault>(inc<_size,__FILE__,__LINE__);
   if (inc!=_inc)
   {
      _inc = inc;
      _base->read(_val.bytes,_ptr + _inc*_val.size,_val.size);
   }
   return _val;
}



template<class M,class T> inline void FileMem::RawPtr<M,T>::operator++()
{
   if (_inc!=(_size - 1))
   {
      _base->read(_val.bytes,_ptr + (++_inc)*_val.size,_val.size);
   }
}



template<class M,class T> inline void FileMem::RawPtr<M,T>::operator--()
{
   if (_inc!=0)
   {
      _base->read(_val.bytes,_ptr + (--_inc)*_val.size,_val.size);
   }
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
