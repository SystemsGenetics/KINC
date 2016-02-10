#ifndef LINUXFILE_H
#define LINUXFILE_H
#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif
#include <string>
#include <cstdint>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "exception.h"



class LinuxFile
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
   // *
   // * BASIC METHODS
   // *
   LinuxFile(const LinuxFile&) = delete;
   LinuxFile(LinuxFile&&) = delete;
   LinuxFile& operator=(const LinuxFile&) = delete;
   LinuxFile& operator=(LinuxFile&&) = delete;
   LinuxFile(const std::string&);
   inline ~LinuxFile();
   // *
   // * FUNCTIONS
   // *
   void clear();
   inline bool reserve(SizeT);
   inline SizeT size() const;
   inline SizeT capacity() const;
   inline Ptr head();
   inline Ptr allocate(SizeT);
   // *
   // * CONSTANTS
   // *
   constexpr static Ptr nullPtr {0xffffffffffffffffull};
   constexpr const static char* _identString = "\0\15\41\102\104\101\124\0\0";
   constexpr static int _idLen = 9;
protected:
   // *
   // * FUNCTIONS
   // *
   inline void write(const void*,Ptr,SizeT);
   inline void read(void*,Ptr,SizeT) const;
private:
   // *
   // * VARIABLES
   // *
   /// File descriptor for file memory object.
   int _fd;
   /// Total size of the file memory object in bytes.
   SizeT _size;
   /// Space available for allocation in bytes.
   SizeT _capacity;
   /// The next location in file memory that will be returned when new space is
   /// allocated.
   Ptr _next;
};



/// @brief Add additional space to object.
///
/// Attempts to add additional space to file memory object that can be used for
/// allocation. If successful at increasing the size of file then additional
/// space is added to available space for allocation.
///
/// @param newBytes Number of additional bytes to add to object.
/// @return True if additional space added.
inline bool LinuxFile::reserve(SizeT size)
{
   bool ret = false;
   if (posix_fallocate64(_fd,lseek64(_fd,0,SEEK_END),size)==0)
   {
      ret = true;
      _size += size;
      _capacity += size;
   }
   return ret;
}



/// If the file descriptor is valid then it is closed.
inline LinuxFile::~LinuxFile()
{
   if (_fd!=-1)
   {
      close(_fd);
   }
}



/// Get total size of file memory object.
///
/// @return Size of object.
inline LinuxFile::SizeT LinuxFile::size() const
{
   return _size;
}



/// Get available space for allocation of file memory object.
///
/// @return Available space of object.
inline LinuxFile::SizeT LinuxFile::capacity() const
{
   return _capacity;
}



/// Get beginning of file memory object.
///
/// @return Points to beginning of file memory.
inline LinuxFile::Ptr LinuxFile::head()
{
   return 0;
}



/// Allocate new space from file memory object.
///
/// @param size Total size in bytes to be allocated.
/// @return Points to beginning of newly allocated file memory.
///
/// @pre The size of the allocation cannot exceed the total space available
/// for allocation in the file memory object.
inline LinuxFile::Ptr LinuxFile::allocate(SizeT size)
{
   assert<OutOfMemory>(size<=_capacity,__FILE__,__LINE__);
   Ptr ret {_next};
   _next += size;
   _capacity -= size;
   bool cond = lseek64(_fd,_idLen,SEEK_SET)==_idLen;
   assert<SystemError>(cond,__FILE__,__LINE__,"lseek64");
   cond = ::write(_fd,&_next,sizeof(Ptr))==sizeof(Ptr);
   assert<SystemError>(cond,__FILE__,__LINE__,"write");
   return ret;
}



/// Write block of binary data to file memory object.
///
/// @param data Binary data in memory to write to file.
/// @param ptr Location in file memory that will be written.
/// @param size Size of binary data to write.
///
/// @pre The block location in file memory cannot exceed the size of the file
/// memory object.
inline void LinuxFile::write(const void* data, Ptr ptr, SizeT size)
{
   assert<FileSegFault>((ptr + size)<=_size,__FILE__,__LINE__);
   int64_t seekr = ptr + _idLen + sizeof(Ptr);
   bool cond = lseek64(_fd,seekr,SEEK_SET)==seekr;
   assert<SystemError>(cond,__FILE__,__LINE__,"lseek64");
   cond = ::write(_fd,data,size)==size;
   assert<SystemError>(cond,__FILE__,__LINE__,"write");
}



/// Read block of binary data from file memory object.
///
/// @param data Memory that will be written.
/// @param ptr Location in file memory that will be read from.
/// @param size Size of binary data to read.
///
/// @pre The block location in file memory cannot exceed the size of the file
/// memory object.
inline void LinuxFile::read(void* data, Ptr ptr, SizeT size) const
{
   assert<FileSegFault>((ptr + size)<=_size,__FILE__,__LINE__);
   int64_t seekr = ptr + _idLen + sizeof(Ptr);
   bool cond = lseek64(_fd,seekr,SEEK_SET)==seekr;
   assert<SystemError>(cond,__FILE__,__LINE__,"lseek64");
   cond = ::read(_fd,data,size)==size;
   assert<SystemError>(cond,__FILE__,__LINE__,"read");
}



/// Generic base exception class for all exceptions thrown in LinuxFile class.
struct LinuxFile::Exception : public ::Exception
{
   using ::Exception::Exception;
};

/// Exception that is thrown when a system error occurs.
struct LinuxFile::SystemError : public ::SystemError
{
   using ::SystemError::SystemError;
};

/// Exception that is thrown when a file that is not a valid file memory object
/// is opened.
struct LinuxFile::InvalidFile : public LinuxFile::Exception
{
   InvalidFile(const char* file, int line):
      Exception(file,line,"FileMem::InvalidFile")
   {}
};

/// Exception that is thrown when a read or write operation is attempted that
/// goes beyond the limits of the total file object's size.
struct LinuxFile::FileSegFault : public LinuxFile::Exception
{
   FileSegFault(const char* file, int line):
      Exception(file,line,"FileMem::FileSegFault")
   {}
};

/// Exception that is thrown when an allocation of new file memory is attempted
/// that is greater than the available amount of bytes that can be allocated.
struct LinuxFile::OutOfMemory : public LinuxFile::Exception
{
   OutOfMemory(const char* file, int line):
      Exception(file,line,"FileMem::OutOfMemory")
   {}
};



#endif
