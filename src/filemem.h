#ifndef FILEMEM_H
#define FILEMEM_H
#include <string>
#include <cstdint>
#include <unistd.h>
#include "exception.h"



enum class FileSync { read,write };



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
   using Ptr = uint64_t;
   using SizeT = uint64_t;
   template<int S> struct Static;
   struct Object;
   class Map;
   // *
   // * BASIC METHODS
   // *
   FileMem(const std::string&);
   ~FileMem();
   // *
   // * FUNCTIONS
   // *
   void clear();
   bool reserve(SizeT);
   inline SizeT size() const;
   inline SizeT capacity() const;
   template<class T> inline void allocate(T&,SizeT = 1);
   template<class T> inline void sync(T&,FileSync,Ptr = 0);
   Ptr head();
   // *
   // * CONSTANTS
   // *
   constexpr static Ptr nullPtr = 0xffffffffffffffffll;
   constexpr static int PtrSize = sizeof(Ptr);
   constexpr const static char* _identString = "\0\15\41\102\104\101\124\0\0";
   constexpr static int _idLen = 9;
private:
   // *
   // * FUNCTIONS
   // *
   inline void write(const void* data,Ptr,SizeT);
   inline void read(void*,Ptr,SizeT) const;
   // *
   // * VARIABLES
   // *
   int _fd;
   SizeT _size;
   SizeT _capacity;
   Ptr _next;
};



template<int S> struct FileMem::Static
{
   // *
   // * DECLERATIONS
   // *
   friend class FileMem;
   // *
   // * BASIC METHODS
   // *
   Static() = default;
   Static(Ptr);
   inline void operator=(Ptr);
   // *
   // * FUNCTIONS
   // *
   template<class T> inline T& get(int);
   inline Ptr addr();
private:
   // *
   // * VARIABLES
   // *
   char bytes[S];
   Ptr ptr {FileMem::nullPtr};
   // *
   // * CONSTANTS
   // *
   constexpr static SizeT size {S};
};



struct FileMem::Object
{
   // *
   // * DECLERATIONS
   // *
   friend class FileMem;
   // *
   // * BASIC METHODS
   // *
   inline Object(const Object&);
   inline Object(Object&&);
   inline Object& operator=(const Object&);
   inline Object& operator=(Object&&);
   inline Object(SizeT,Ptr = nullPtr);
   inline ~Object();
   inline void operator=(Ptr);
   // *
   // * FUNCTIONS
   // *
   template<class T> inline T& get(int);
   inline Ptr addr();
private:
   // *
   // * VARIABLES
   // *
   char* bytes;
   SizeT size;
   Ptr ptr {FileMem::nullPtr};
};



inline FileMem::SizeT FileMem::size() const
{
   return _size;
}



inline FileMem::SizeT FileMem::capacity() const
{
   return _capacity;
}



template<class T> inline void FileMem::allocate(T& obj, SizeT amt)
{
   if (amt>0)
   {
      SizeT size = obj.size*amt;
      assert<OutOfMemory>(size<=_capacity,__FILE__,__LINE__);
      obj.ptr = _next;
      _next += size;
      _capacity -= size;
      bool cond = lseek64(_fd,_idLen,SEEK_SET)==_idLen;
      assert<SystemError>(cond,__FILE__,__LINE__,"lseek64");
      cond = ::write(_fd,&_next,sizeof(Ptr))==sizeof(Ptr);
      assert<SystemError>(cond,__FILE__,__LINE__,"write");
   }
}



template<class T> inline void FileMem::sync(T& obj, FileSync which, Ptr inc)
{
   assert<FileSegFault>(obj.ptr!=nullPtr,__FILE__,__LINE__);
   Ptr seekr = obj.ptr + obj.size*inc;
   switch (which)
   {
   case FileSync::read:
      read(obj.bytes,seekr,obj.size);
      break;
   case FileSync::write:
      write(obj.bytes,seekr,obj.size);
      break;
   }
}



inline FileMem::Ptr FileMem::head()
{
   return 0;
}



inline void FileMem::write(const void* data, Ptr ptr, SizeT size)
{
   assert<FileSegFault>((ptr + size)<=_size,__FILE__,__LINE__);
   int64_t seekr = ptr + _idLen + sizeof(Ptr);
   bool cond = lseek64(_fd,seekr,SEEK_SET)==seekr;
   assert<SystemError>(cond,__FILE__,__LINE__,"lseek64");
   cond = ::write(_fd,data,size)==size;
   assert<SystemError>(cond,__FILE__,__LINE__,"write");
}



inline void FileMem::read(void* data, Ptr ptr, SizeT size) const
{
   assert<FileSegFault>((ptr + size)<=_size,__FILE__,__LINE__);
   int64_t seekr = ptr + _idLen + sizeof(Ptr);
   bool cond = lseek64(_fd,seekr,SEEK_SET)==seekr;
   assert<SystemError>(cond,__FILE__,__LINE__,"lseek64");
   cond = ::read(_fd,data,size)==size;
   assert<SystemError>(cond,__FILE__,__LINE__,"read");
}



template<int S> FileMem::Static<S>::Static(Ptr p):
   ptr(p)
{}



template<int S> inline void FileMem::Static<S>::operator=(Ptr p)
{
   ptr = p;
}



template<int S> template<class T> inline T& FileMem::Static<S>::get(int n)
{
   return *reinterpret_cast<T*>(&bytes[n]);
}



template<int S> inline FileMem::Ptr FileMem::Static<S>::addr()
{
   return ptr;
}



inline FileMem::Object::Object(const Object& obj):
   size(obj.size)
{
   bytes = new char[size];
   memcpy(bytes,obj.bytes,size);
}



inline FileMem::Object::Object(Object&& tmp):
   size(tmp.size),
   bytes(tmp.bytes)
{
   tmp.bytes = nullptr;
   tmp.size = 0;
}



inline FileMem::Object& FileMem::Object::operator=(const Object& obj)
{
   size = obj.size;
   if (bytes)
   {
      delete[] bytes;
   }
   bytes = new char[size];
   memcpy(bytes,obj.bytes,size);
}



inline FileMem::Object& FileMem::Object::operator=(Object&& tmp)
{
   size = tmp.size;
   bytes = tmp.bytes;
   tmp.bytes = nullptr;
   tmp.size = 0;
}



inline FileMem::Object::Object(SizeT s, Ptr p):
   size(s),
   ptr(p)
{
   bytes = new char[size];
}



inline FileMem::Object::~Object()
{
   if (bytes)
   {
      delete[] bytes;
   }
}



inline void FileMem::Object::operator=(Ptr p)
{
   ptr = p;
}



template<class T> inline T& FileMem::Object::get(int n)
{
   return *reinterpret_cast<T*>(&bytes[n]);
}



inline FileMem::Ptr FileMem::Object::addr()
{
   return ptr;
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
