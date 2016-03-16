#include <sys/types.h>
#include <sys/stat.h>
#include <cstring>
#include "filemem.h"
#include <iostream>



FileMem::FileMem(const std::string& fileName):
   _fd(-1),
   _size(0),
   _capacity(0),
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
             ::write(_fd,&_next,sizeof(Ptr))==sizeof(Ptr);
      assert<SystemError>(cond,__FILE__,__LINE__,"write");
   }
   else
   {
      cond = lseek64(_fd,0,SEEK_END)>=(_idLen + sizeof(Ptr));
      assert<InvalidFile>(cond,__FILE__,__LINE__);
      char buf[_idLen];
      cond = ::lseek64(_fd,0,SEEK_SET)==0;
      assert<SystemError>(cond,__FILE__,__LINE__,"lseek64");
      cond = ::read(_fd,buf,_idLen)!=-1;
      assert<SystemError>(cond,__FILE__,__LINE__,"read");
      cond = true;
      for (int i=0;i<_idLen;++i)
      {
         if (buf[i]!=_identString[i])
         {
            cond = false;
         }
      }
      assert<InvalidFile>(cond,__FILE__,__LINE__);
      cond = ::read(_fd,&_next,sizeof(Ptr))==sizeof(Ptr);
      assert<SystemError>(cond,__FILE__,__LINE__,"read");
      _size = lseek64(_fd,0,SEEK_END) - _idLen - sizeof(Ptr);
      _capacity = _size - _next;
   }
}



FileMem::~FileMem()
{
   if (_fd!=-1)
   {
      close(_fd);
   }
}



void FileMem::clear()
{
   _capacity = _size;
   _next = 0;
   bool cond = lseek64(_fd,_idLen,SEEK_SET)==_idLen;
   assert<SystemError>(cond,__FILE__,__LINE__,"lseek64");
   cond = ::write(_fd,&_next,sizeof(Ptr))==sizeof(Ptr);
   assert<SystemError>(cond,__FILE__,__LINE__,"write");
}
