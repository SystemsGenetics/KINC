#include <utility>
#include "filemem.h"
#include "unit.h"



namespace unit
{
   namespace filemem
   {
      struct Node : FileMem::Static<12>
      {
         using Static<12>::Static;
         Node()
         {
            val() = 42;
         }
         uint32_t& val() { get<uint32_t>(0); }
         FileMem::Ptr& next() { get<FileMem::Ptr>(4); }
      };
      struct LNode : FileMem::Static<5285>
      {};
      struct DNode : FileMem::Object
      {
         using Object::Object;
         DNode(): Object(12)
         {
            val() = 42;
         }
         DNode(FileMem::Ptr p): Object(12,p)
         {
            val() = 42;
         }
         uint32_t& val() { get<uint32_t>(0); }
         FileMem::Ptr& next() { get<FileMem::Ptr>(4); }
      };
      constexpr const char* tmpFile = "memfile.tmp";
      constexpr const char* validFile = "memfile2.tmp";
      constexpr const char* invalidFile = "notmemfile.tmp";
      constexpr const char* invalidFile2 = "notmemfile2.tmp";
      constexpr const char* wrtData = "1234567890123456";
   }
}



void construct_filemems()
{
   constexpr static int flags = O_CREAT|O_RDWR|O_LARGEFILE;
   constexpr static mode_t modes = S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH;
   int fd;
   char identString[10] = "\0\15\41\102\104\101\124\0\0";
   int idLen = 9;
   FileMem::Ptr next = 0;
   fd = open(unit::filemem::validFile,flags,modes);
   ::write(fd,identString,idLen);
   ::write(fd,&next,sizeof(FileMem::Ptr));
   close(fd);
   fd = open(unit::filemem::invalidFile,flags,modes);
   ::write(fd,identString,idLen);
   ::write(fd,&next,sizeof(FileMem::Ptr)-1);
   close(fd);
   identString[3] = '\16';
   fd = open(unit::filemem::invalidFile2,flags,modes);
   ::write(fd,identString,idLen);
   ::write(fd,&next,sizeof(FileMem::Ptr));
   close(fd);
}



bool unit::filemem::main()
{
   bool ret = false;
   header("FileMem");
   try
   {
      construct_filemems();
      ret = static1()&&
            object1()&&
            object2()&&
            object3()&&
            object4()&&
            object5()&&
            init1()&&
            init2()&&
            init3()&&
            init4()&&
            reserve1()&&
            expand1()&&
            capacity1()&&
            allocate1()&&
            allocate2()&&
            allocate3()&&
            allot1()&&
            clear1()&&
            addr1()&&
            sync1()&&
            sync2()&&
            sync3()&&
            sync4()&&
            sync5()&&
            sync6();
   }
   catch (...)
   {
      system("rm -f *.tmp");
      end();
      throw;
   }
   system("rm -f *.tmp");
   end();
   return ret;
}



bool unit::filemem::static1()
{
   start();
   Node t;
   bool ret = t.val()==42;
   return finish(ret,"static1");
}



bool unit::filemem::object1()
{
   start();
   DNode t;
   bool ret = t.val()==42;
   return finish(ret,"object1");
}



bool unit::filemem::object2()
{
   start();
   DNode cpy;
   DNode t(cpy);
   bool ret = t.val()==42;
   return finish(ret,"object2");
}



bool unit::filemem::object3()
{
   start();
   DNode cpy;
   DNode t(std::move(cpy));
   bool ret = t.val()==42;
   return finish(ret,"object3");
}



bool unit::filemem::object4()
{
   start();
   DNode cpy;
   DNode t;
   t = cpy;
   bool ret = t.val()==42;
   return finish(ret,"object4");
}



bool unit::filemem::object5()
{
   start();
   DNode cpy;
   DNode t;
   t = std::move(cpy);
   bool ret = t.val()==42;
   return finish(ret,"object5");
}



bool unit::filemem::init1()
{
   start();
   FileMem tf(tmpFile);
   bool ret = tf.size()==0;
   return finish(ret,"init1");
}



bool unit::filemem::init2()
{
   start();
   FileMem tf(validFile);
   bool ret = tf.size()==0;
   return finish(ret,"init2");
}



bool unit::filemem::init3()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(invalidFile);
      tf.size();
   }
   catch (FileMem::InvalidFile)
   {
      ret = true;
   }
   return finish(ret,"init3");
}



bool unit::filemem::init4()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(invalidFile2);
      tf.size();
   }
   catch (FileMem::InvalidFile)
   {
      ret = true;
   }
   return finish(ret,"init4");
}



bool unit::filemem::reserve1()
{
   start();
   FileMem tf(tmpFile);
   tf.reserve(4096);
   bool ret = tf.size()==4096;
   return finish(ret,"reserve1");
}



bool unit::filemem::expand1()
{
   start();
   FileMem tf(tmpFile);
   Node tmp;
   tf.expand(tmp,100);
   bool ret = tf.size()==5296;
   return finish(ret,"expand1");
}



bool unit::filemem::capacity1()
{
   start();
   FileMem tf(tmpFile);
   bool ret = tf.capacity()==5296;
   return finish(ret,"capacity1");
}



bool unit::filemem::allocate1()
{
   start();
   FileMem tf(tmpFile);
   Node tmp;
   tf.allocate(tmp);
   bool ret = tf.capacity()==5284;
   return finish(ret,"allocate1");
}



bool unit::filemem::allocate2()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      LNode tmp;
      tf.allocate(tmp);
   }
   catch (FileMem::OutOfMemory)
   {
      ret = true;
   }
   return finish(ret,"allocate2");
}



bool unit::filemem::allocate3()
{
   start();
   bool ret = true;
   try
   {
      FileMem tf(tmpFile);
      LNode tmp;
      tf.allocate(tmp,0);
   }
   catch (FileMem::OutOfMemory)
   {
      ret = false;
   }
   return finish(ret,"allocate3");
}



bool unit::filemem::allot1()
{
   start();
   bool ret = true;
   try
   {
      FileMem tf(tmpFile);
      LNode tmp;
      tf.allot(tmp);
   }
   catch (FileMem::OutOfMemory)
   {
      ret = false;
   }
   return finish(ret,"allot1");
}



bool unit::filemem::clear1()
{
   start();
   FileMem tf(tmpFile);
   tf.clear();
   Node tmp;
   tf.allocate(tmp);
   bool ret {tmp.addr()==tf.head()};
   return finish(ret,"clear1");
}



bool unit::filemem::addr1()
{
   start();
   FileMem tf(tmpFile);
   tf.clear();
   DNode tmp;
   tf.allocate(tmp);
   bool ret {tmp.addr()==tf.head()};
   return finish(ret,"addr1");
}



bool unit::filemem::sync1()
{
   start();
   {
      FileMem tf(tmpFile);
      tf.clear();
      Node data;
      data.val() = 33;
      tf.allocate(data);
      tf.sync(data,FileSync::write);
   }
   FileMem tf(tmpFile);
   Node data = tf.head();
   tf.sync(data,FileSync::read);
   bool ret {data.val()==33};
   return finish(ret,"sync1");
}



bool unit::filemem::sync2()
{
   start();
   {
      FileMem tf(tmpFile);
      tf.clear();
      Node data;
      data.val() = 33;
      tf.allocate(data,3);
      tf.sync(data,FileSync::write,2);
   }
   FileMem tf(tmpFile);
   Node data = tf.head();
   tf.sync(data,FileSync::read,2);
   bool ret {data.val()==33};
   return finish(ret,"sync2");
}



bool unit::filemem::sync3()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      Node data;
      tf.sync(data,FileSync::read);
   }
   catch (FileMem::FileSegFault)
   {
      ret = true;
   }
   return finish(ret,"sync3");
}



bool unit::filemem::sync4()
{
   start();
   {
      FileMem tf(tmpFile);
      tf.clear();
      DNode data;
      data.val() = 33;
      tf.allocate(data);
      tf.sync(data,FileSync::write);
   }
   FileMem tf(tmpFile);
   Node data = tf.head();
   tf.sync(data,FileSync::read);
   bool ret {data.val()==33};
   return finish(ret,"sync4");
}



bool unit::filemem::sync5()
{
   start();
   {
      FileMem tf(tmpFile);
      tf.clear();
      DNode data;
      data.val() = 33;
      tf.allocate(data,3);
      tf.sync(data,FileSync::write,2);
   }
   FileMem tf(tmpFile);
   DNode data = tf.head();
   tf.sync(data,FileSync::read,2);
   bool ret {data.val()==33};
   return finish(ret,"sync5");
}



bool unit::filemem::sync6()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      DNode data;
      tf.sync(data,FileSync::read);
   }
   catch (FileMem::FileSegFault)
   {
      ret = true;
   }
   return finish(ret,"sync6");
}
