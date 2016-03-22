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
   char identString[6] = "\33\102\104\101\124";
   int idLen = 5;
   FileMem::Ptr next = 0;
   fd = open(unit::filemem::validFile,flags,modes);
   ::write(fd,identString,idLen);
   ::write(fd,&next,sizeof(FileMem::Ptr));
   close(fd);
   fd = open(unit::filemem::invalidFile,flags,modes);
   ::write(fd,identString,idLen);
   ::write(fd,&next,sizeof(FileMem::Ptr)-1);
   close(fd);
   identString[1] = '\16';
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
      ret = static_()&&
            object()&&
            construct()&&
            reserve()&&
            expand()&&
            allocate()&&
            allot()&&
            clear()&&
            addr()&&
            sync();
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



bool unit::filemem::static_()
{
   start();
   Node t;
   bool test = t.val()==42;
   return finish (test,"static1");
}



bool unit::filemem::object()
{
   bool cont = true;
   {
      start();
      DNode t;
      bool test = t.val()==42;
      cont = cont&&finish (test,"object1");
   }
   if (cont)
   {
      start();
      DNode cpy;
      DNode t(cpy);
      bool test = t.val()==42;
      cont = cont&&finish (test,"object2");
   }
   if (cont)
   {
      start();
      DNode cpy;
      DNode t(std::move(cpy));
      bool test = t.val()==42;
      cont = cont&&finish (test,"object3");
   }
   if (cont)
   {
      start();
      DNode cpy;
      DNode t;
      t = cpy;
      bool test = t.val()==42;
      cont = cont&&finish (test,"object4");
   }
   if (cont)
   {
      start();
      DNode cpy;
      DNode t;
      t = std::move(cpy);
      bool test = t.val()==42;
      cont = cont&&finish (test,"object5");
   }
   return cont;
}



bool unit::filemem::construct()
{
   bool cont = true;
   {
      start();
      FileMem tf(tmpFile);
      bool test = tf.size()==0;
      cont = cont&&finish (test,"construct1");
   }
   {
      start();
      FileMem tf(validFile);
      bool test = tf.size()==0;
      cont = cont&&finish (test,"construct2");
   }
   {
      start();
      bool test = false;
      try
      {
         FileMem tf(invalidFile);
      }
      catch (FileMem::InvalidFile)
      {
         test = true;
      }
      cont = cont&&finish (test,"construct3");
   }
   {
      start();
      bool test = false;
      try
      {
         FileMem tf(invalidFile2);
      }
      catch (FileMem::InvalidFile)
      {
         test = true;
      }
      cont = cont&&finish (test,"construct4");
   }
   return cont;
}



bool unit::filemem::reserve()
{
   start();
   FileMem tf(tmpFile);
   tf.reserve(4096);
   bool test = tf.capacity()==4096;
   return finish (test,"reserve1");
}



bool unit::filemem::expand()
{
   start();
   FileMem tf(tmpFile);
   Node tmp;
   tf.expand(tmp,100);
   bool test = tf.capacity()==5296;
   return finish (test,"expand1");
}



bool unit::filemem::allocate()
{
   bool cont = true;
   {
      start();
      FileMem tf(tmpFile);
      Node tmp;
      tf.allocate(tmp);
      bool test = tf.available()==5284;
      cont = cont&&finish (test,"allocate1");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      Node tmp;
      tf.allocate(tmp);
      bool test = tf.size()==24;
      cont = cont&&finish (test,"allocate2");
   }
   if (cont)
   {
      start();
      bool test = false;
      try
      {
         FileMem tf(tmpFile);
         LNode tmp;
         tf.allocate(tmp);
      }
      catch (FileMem::OutOfMemory)
      {
         test = true;
      }
      cont = cont&&finish (test,"allocate3");
   }
   if (cont)
   {
      start();
      bool test = true;
      try
      {
         FileMem tf(tmpFile);
         LNode tmp;
         tf.allocate(tmp,0);
      }
      catch (FileMem::OutOfMemory)
      {
         test = false;
      }
      cont = cont&&finish (test,"allocate4");
   }
   return cont;
}



bool unit::filemem::allot()
{
   start();
   bool test = true;
   try
   {
      FileMem tf(tmpFile);
      LNode tmp;
      tf.allot(tmp);
   }
   catch (FileMem::OutOfMemory)
   {
      test = false;
   }
   return finish (test,"allot1");
}



bool unit::filemem::clear()
{
   start();
   FileMem tf(tmpFile);
   tf.clear();
   Node tmp;
   tf.allocate(tmp);
   bool test {tmp.addr()==tf.head()};
   return finish (test,"clear1");
}



bool unit::filemem::addr()
{
   start();
   FileMem tf(tmpFile);
   tf.clear();
   DNode tmp;
   tf.allocate(tmp);
   bool test {tmp.addr()==tf.head()};
   return finish (test,"addr1");
}



bool unit::filemem::sync()
{
   bool cont {true};
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
      bool test {data.val()==33};
      cont = cont&&finish (test,"sync1");
   }
   if (cont)
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
      bool test {data.val()==33};
      cont = cont&&finish (test,"sync2");
   }
   if (cont)
   {
      start();
      bool test = false;
      try
      {
         FileMem tf(tmpFile);
         Node data;
         tf.sync(data,FileSync::read);
      }
      catch (FileMem::FileSegFault)
      {
         test = true;
      }
      cont = cont&&finish (test,"sync3");
   }
   if (cont)
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
      bool test {data.val()==33};
      cont = cont&&finish (test,"sync4");
   }
   if (cont)
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
      bool test {data.val()==33};
      cont = cont&&finish (test,"sync5");
   }
   if (cont)
   {
      start();
      bool test = false;
      try
      {
         FileMem tf(tmpFile);
         DNode data;
         tf.sync(data,FileSync::read);
      }
      catch (FileMem::FileSegFault)
      {
         test = true;
      }
      cont = cont&&finish (test,"sync6");
   }
   if (cont)
   {
      start();
      bool test = false;
      try
      {
         FileMem tf(tmpFile);
         tf.clear();
         DNode data(tf.head());
         tf.sync(data,FileSync::read);
      }
      catch (FileMem::FileSegFault)
      {
         test = true;
      }
      cont = cont&&finish (test,"sync7");
   }
   if (cont)
   {
      start();
      bool test = false;
      try
      {
         FileMem tf(tmpFile);
         DNode data(tf.head());
         tf.sync(data,FileSync::write);
      }
      catch (FileMem::FileSegFault)
      {
         test = true;
      }
      cont = cont&&finish (test,"sync8");
   }
   return cont;
}
