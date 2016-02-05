#include "linuxfile.h"



/// @brief Opens memory object file and initializes object.
///
/// Opens the memory object file. If the file is new(0 size) then it writes
/// to the file making a new memory object with 0 size. If it is not new(size
/// greater than 0) then it reads the header of the file to make sure it is a
/// valid memory object file. If it is valid, it reads in the next pointer and
/// calculates total size and available bytes of file object.
///
/// @pre The file opened must be a valid memory object or a new file.
LinuxFile::LinuxFile(const std::string& fileName):
   _fd(-1),
   _size(0),
   _available(0),
   _next(0)
{
   constexpr static int flags = O_CREAT|O_RDWR|O_LARGEFILE;
   constexpr static mode_t modes = S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH;
   _fd = open(fileName.c_str(),flags,modes);
   assert<SystemError>(_fd!=-1,__FILE__,__LINE__,"open");
   bool cond;
   if (lseek64(_fd,0,SEEK_END)==0)
   {
      cond = ::write(_fd,_identString,_idLen)==_idLen&&
             ::write(_fd,&_next,sizeof(VPtr))==sizeof(VPtr);
      assert<SystemError>(cond,__FILE__,__LINE__,"write");
   }
   else
   {
      char buf[_idLen];
      cond = ::lseek64(_fd,0,SEEK_SET)==0;
      assert<SystemError>(cond,__FILE__,__LINE__,"lseek64");
      cond = ::read(_fd,buf,_idLen)!=-1;
      assert<SystemError>(cond,__FILE__,__LINE__,"read");
      cond = true;
      for (int i=0;i<_idLen;i++)
      {
         if (buf[i]!=_identString[i])
         {
            cond = false;
         }
      }
      assert<InvalidFile>(cond,__FILE__,__LINE__);
      cond = ::read(_fd,&_next,sizeof(VPtr))==sizeof(VPtr);
      assert<SystemError>(cond,__FILE__,__LINE__,"read");
      _size = lseek64(_fd,0,SEEK_END) - _idLen - sizeof(VPtr);
      _available = _size - _next;
   }
}



/// @brief Clears allocated memory of object.
///
/// Clears the file memory object of all allocated space, returning the total
/// size of the object as available bytes that can be allocated.
void LinuxFile::clear()
{
   _available = _size;
   _next = 0;
   bool cond = lseek64(_fd,_idLen,SEEK_SET)==_idLen;
   assert<SystemError>(cond,__FILE__,__LINE__,"lseek64");
   cond = ::write(_fd,&_next,sizeof(VPtr))==sizeof(VPtr);
   assert<SystemError>(cond,__FILE__,__LINE__,"write");
}



/// @brief Add additional space to object.
///
/// Attempts to add additional space to file memory object that can be used for
/// allocation. If successful at increasing the size of file then additional
/// space is added to available space for allocation.
///
/// @param newBytes Number of additional bytes to add to object.
/// @return True if additional space added.
bool LinuxFile::reserve(int64_t newBytes)
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
