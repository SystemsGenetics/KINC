#ifndef FILEMEM_H
#define FILEMEM_H
#include <string>
#include <cstdint>
#include "exception.h"



class FileMem
{
public:
   struct Exception;
   struct SystemError;
   struct InvalidFile;
   struct FileSegFault;
   struct OutOfMemory;
   template<class T> class Ptr;
   using VPtr = uint64_t;
   virtual ~FileMem() = default;
   virtual void clear() = 0;
   virtual bool reserve(int64_t) = 0;
   virtual uint64_t size() = 0;
   virtual uint64_t available() = 0;
   constexpr const static char* _identString = "KINCBINDAT";
protected:
   virtual void write(const void*,VPtr,uint64_t) = 0;
   virtual void read(void*,VPtr,uint64_t) = 0;
   virtual VPtr allocate(uint64_t) = 0;
   virtual VPtr head() = 0;
};



// in <type_traits> std::is_integral<T>::value &&
//                  std::is_floating_point<T>::value :D
// and std::is_class<T>::value AND std::is_polymorphic<T>::value :D
template<class T> class FileMem::Ptr
{
public:
   Ptr(const Ptr<T>&) = delete;
   Ptr<T>& operator=(const Ptr<T>&) = delete;
   Ptr(FileMem&,T,int);
   Ptr(FileMem&,int);
   Ptr(FileMem&,VPtr,int);
   Ptr(Ptr<T>&&);
   Ptr<T>& operator=(Ptr<T>&&);
   Ptr<T>& operator=(VPtr);
   T& operator*();
   T& operator[](int);
   VPtr addr();
   void save();
   void free();
};



struct FileMem::Exception : public ::Exception
{
   using ::Exception::Exception;
};

struct FileMem::SystemError : public ::SystemError
{
   using ::SystemError::SystemError;
};

struct FileMem::InvalidFile : public FileMem::Exception
{
   InvalidFile(const char* file, int line):
      Exception(file,line,"FileMem::InvalidFile")
   {}
};

struct FileMem::FileSegFault : public FileMem::Exception
{
   FileSegFault(const char* file, int line):
      Exception(file,line,"FileMem::FileSegFault")
   {}
};

struct FileMem::OutOfMemory : public FileMem::Exception
{
   OutOfMemory(const char* file, int line):
      Exception(file,line,"FileMem::OutOfMemory")
   {}
};



#endif
