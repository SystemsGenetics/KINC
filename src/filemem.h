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
   char bytes[S];
   constexpr static uint64_t size = S;
   template<class T> T& get(int n)
   {
      return *reinterpret_cast<T*>(&bytes[n]);
   }
};



struct FileMem::DObject
{
   DObject(uint64_t dSize):
      size(dSize)
   {
      bytes = new char[size];
   }
   DObject(const DObject& obj):
      size(obj.size)
   {
      bytes = new char[size];
      memcpy(bytes,obj.bytes,size);
   }
   DObject(DObject&& tmp):
      size(tmp.size),
      bytes(tmp.bytes)
   {
      tmp.bytes = nullptr;
      tmp.size = 0;
   }
   DObject& operator=(const DObject& obj)
   {
      size = obj.size;
      if (bytes)
      {
         delete[] bytes;
      }
      bytes = new char[size];
      memcpy(bytes,obj.bytes,size);
   }
   DObject& operator=(DObject&& tmp)
   {
      size = tmp.size;
      bytes = tmp.bytes;
      tmp.bytes = nullptr;
      tmp.size = 0;
   }
   ~DObject()
   {
      if (bytes)
      {
         delete[] bytes;
      }
   }
   char* bytes;
   uint64_t size;
   template<class T> T& get(int n)
   {
      return *reinterpret_cast<T*>(&bytes[n]);
   }
};



template<class M,class T> class FileMem::RawPtr
{
public:
   static_assert(std::is_class<T>::value,
                 "FileMem::Ptr<T>, T is of an unsupported type.");
   RawPtr(const Ptr<T>&) = delete;
   RawPtr<M,T>& operator=(const Ptr<T>&) = delete;
   RawPtr(M&);
   RawPtr(M&,const T&);
   RawPtr(M&,T&&);
   RawPtr(M&,VPtr);
   RawPtr(RawPtr<M,T>&&);
   inline void operator=(RawPtr<M,T>&&);
   inline void operator=(VPtr);
   inline T& operator*();
   inline const T& operator*() const;
   inline T* operator->();
   inline const T* operator->() const;
   inline VPtr addr() const;
   inline void save();
private:
   T _val;
   M* _base;
   VPtr _ptr;
};



template<class M,class T> FileMem::RawPtr<M,T>::RawPtr(M& mem):
   _base(&mem),
   _ptr(nullPtr)
{
   _ptr = _base->allocate(_val.size);
}



template<class M,class T> FileMem::RawPtr<M,T>::RawPtr(M& mem, const T& val):
   _val(val),
   _base(&mem),
   _ptr(nullPtr)
{
   _ptr = _base->allocate(val.size);
}



template<class M,class T> FileMem::RawPtr<M,T>::RawPtr(M& mem, T&& val):
   _val(std::move(val)),
   _base(&mem),
   _ptr(nullPtr)
{
   _ptr = _base->allocate(val.size);
}



template<class M,class T> FileMem::RawPtr<M,T>::RawPtr(M& mem, VPtr loc):
   _base(&mem),
   _ptr(loc)
{
   _base->read(_val.bytes,_ptr,_val.size);
}



template<class M,class T>
FileMem::RawPtr<M,T>::RawPtr(FileMem::RawPtr<M,T>&& tmp):
   _val(std::move(tmp._val)),
   _base(tmp._base),
   _ptr(tmp._ptr)
{
   tmp._ptr = nullPtr;
}



template<class M,class T>
inline void FileMem::RawPtr<M,T>::operator=(RawPtr<M,T>&& tmp)
{
   _val = std::move(tmp._val);
   _base = tmp._base;
   _ptr = tmp._ptr;
   tmp._ptr = nullPtr;
}



template<class M,class T>
inline void FileMem::RawPtr<M,T>::operator=(VPtr ptr)
{
   _ptr = ptr;
   _base->read(_val.bytes,_ptr,_val.size);
}



template<class M,class T> inline T& FileMem::RawPtr<M,T>::operator*()
{
   return _val;
}



template<class M,class T>
inline const T& FileMem::RawPtr<M,T>::operator*() const
{
   return _val;
}



template<class M,class T> inline T* FileMem::RawPtr<M,T>::operator->()
{
   return &_val;
}



template<class M,class T>
inline const T* FileMem::RawPtr<M,T>::operator->() const
{
   return &_val;
}



template<class M,class T>
inline FileMem::VPtr FileMem::RawPtr<M,T>::addr() const
{
   return _ptr;
}



template<class M,class T> inline void FileMem::RawPtr<M,T>::save()
{
   _base->write(_val.bytes,_ptr,_val.size);
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
