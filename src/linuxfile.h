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
#include "filemem.h"



/// @brief Linux OS FileMem implementation.
///
/// Implements the file memory base class within the Linux OS environment.
/// Uses basic file IO calls in linux such as open/close/read/write and
/// posix_fallocate.
///
/// @author Josh Burns
/// @date 3 Feb 2016
class LinuxFile : public FileMem
{
public:
   // *
   // * DECLERATIONS
   // *
   template<class M,class T> friend class FileMem::RawPtr;
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
   inline bool reserve(int64_t);
   inline uint64_t size() const;
   inline uint64_t available() const;
   inline VPtr head() const;
protected:
   // *
   // * FUNCTIONS
   // *
   inline void write(const void*,VPtr,uint64_t);
   inline void read(void*,VPtr,uint64_t) const;
   inline VPtr allocate(uint64_t);
private:
   // *
   // * VARIABLES
   // *
   /// File descriptor for file memory object.
   int _fd;
   /// Total size of the file memory object in bytes.
   uint64_t _size;
   /// Space available for allocation in bytes.
   uint64_t _available;
   /// The next location in file memory that will be returned when new space is
   /// allocated.
   VPtr _next;
};



/// @brief Add additional space to object.
///
/// Attempts to add additional space to file memory object that can be used for
/// allocation. If successful at increasing the size of file then additional
/// space is added to available space for allocation.
///
/// @param newBytes Number of additional bytes to add to object.
/// @return True if additional space added.
inline bool LinuxFile::reserve(int64_t newBytes)
{
   bool ret = false;
   if (posix_fallocate64(_fd,lseek64(_fd,0,SEEK_END),newBytes)==0)
   {
      ret = true;
      _size += newBytes;
      _available += newBytes;
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
inline uint64_t LinuxFile::size() const
{
   return _size;
}



/// Get available space for allocation of file memory object.
///
/// @return Available space of object.
inline uint64_t LinuxFile::available() const
{
   return _available;
}



/// Get beginning of file memory object.
///
/// @return Points to beginning of file memory.
inline FileMem::VPtr LinuxFile::head() const
{
   return 0;
}



/// Write block of binary data to file memory object.
///
/// @param data Binary data in memory to write to file.
/// @param ptr Location in file memory that will be written.
/// @param size Size of binary data to write.
///
/// @pre The block location in file memory cannot exceed the size of the file
/// memory object.
inline void LinuxFile::write(const void* data, VPtr ptr, uint64_t size)
{
   assert<FileSegFault>((ptr + size)<=_size,__FILE__,__LINE__);
   int64_t seekr = ptr + _idLen + sizeof(VPtr);
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
inline void LinuxFile::read(void* data, VPtr ptr, uint64_t size) const
{
   assert<FileSegFault>((ptr + size)<=_size,__FILE__,__LINE__);
   int64_t seekr = ptr + _idLen + sizeof(VPtr);
   bool cond = lseek64(_fd,seekr,SEEK_SET)==seekr;
   assert<SystemError>(cond,__FILE__,__LINE__,"lseek64");
   cond = ::read(_fd,data,size)==size;
   assert<SystemError>(cond,__FILE__,__LINE__,"read");
}



/// Allocate new space from file memory object.
///
/// @param size Total size in bytes to be allocated.
/// @return Points to beginning of newly allocated file memory.
///
/// @pre The size of the allocation cannot exceed the total space available
/// for allocation in the file memory object.
inline FileMem::VPtr LinuxFile::allocate(uint64_t size)
{
   assert<OutOfMemory>(size<=_available,__FILE__,__LINE__);
   VPtr ret = _next;
   _next += size;
   _available -= size;
   bool cond = lseek64(_fd,_idLen,SEEK_SET)==_idLen;
   assert<SystemError>(cond,__FILE__,__LINE__,"lseek64");
   cond = ::write(_fd,&_next,sizeof(VPtr))==sizeof(VPtr);
   assert<SystemError>(cond,__FILE__,__LINE__,"write");
   return ret;
}



#endif
