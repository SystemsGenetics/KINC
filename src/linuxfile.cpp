#define _LARGEFILE64_SOURCE
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <cstring>
#include <iostream>
#include "linuxfile.h"



LinuxFile::LinuxFile(const std::string& fileName):
   _fd(-1),
   _idLen(strlen(_identString)),
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
      char buf[_idLen+1] {'\0'};
      cond = ::lseek64(_fd,0,SEEK_SET)==0;
      assert<SystemError>(cond,__FILE__,__LINE__,"lseek64");
      cond = ::read(_fd,buf,_idLen)!=-1;
      assert<SystemError>(cond,__FILE__,__LINE__,"read");
      cond = strcmp(buf,_identString)==0;
      assert<InvalidFile>(cond,__FILE__,__LINE__);
      cond = ::read(_fd,&_next,sizeof(VPtr))==sizeof(VPtr);
      assert<SystemError>(cond,__FILE__,__LINE__,"read");
      _size = lseek64(_fd,0,SEEK_END) - _idLen - sizeof(VPtr);
      _available = _size - _next;
   }
}



LinuxFile::~LinuxFile()
{
   if (_fd!=-1)
   {
      close(_fd);
   }
}



void LinuxFile::clear()
{
   _available = _size;
   _next = 0;
   bool cond = lseek64(_fd,_idLen,SEEK_SET)==_idLen;
   assert<SystemError>(cond,__FILE__,__LINE__,"lseek64");
   cond = ::write(_fd,&_next,sizeof(VPtr))==sizeof(VPtr);
   assert<SystemError>(cond,__FILE__,__LINE__,"write");
}



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



uint64_t LinuxFile::size()
{
   return _size;
}



uint64_t LinuxFile::available()
{
   return _available;
}



void LinuxFile::write(const void* data, VPtr ptr, uint64_t size)
{
   assert<FileSegFault>((ptr + size)<=_size,__FILE__,__LINE__);
   int64_t seekr = ptr + _idLen + sizeof(VPtr);
   bool cond = lseek64(_fd,seekr,SEEK_SET)==seekr;
   assert<SystemError>(cond,__FILE__,__LINE__,"lseek64");
   cond = ::write(_fd,data,size)==size;
   assert<SystemError>(cond,__FILE__,__LINE__,"write");
}



void LinuxFile::read(void* data, VPtr ptr, uint64_t size)
{
   assert<FileSegFault>((ptr + size)<=_size,__FILE__,__LINE__);
   int64_t seekr = ptr + _idLen + sizeof(VPtr);
   bool cond = lseek64(_fd,seekr,SEEK_SET)==seekr;
   assert<SystemError>(cond,__FILE__,__LINE__,"lseek64");
   cond = ::read(_fd,data,size)==size;
   assert<SystemError>(cond,__FILE__,__LINE__,"read");
}



FileMem::VPtr LinuxFile::allocate(uint64_t size)
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



FileMem::VPtr LinuxFile::head()
{
   return 0;
}
