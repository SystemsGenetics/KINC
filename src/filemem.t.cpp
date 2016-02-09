#include <cstdlib>
#include <utility>
#include "filemem.h"
#include "unit.h"



struct Node : FileMem::Object<12>
{
   Node()
   {
      val() = 42;
      next() = FileMem::nullPtr;
   }
   uint32_t& val() { get<uint32_t>(0); }
   FileMem::VPtr& next() { get<FileMem::VPtr>(4); }
};



struct DNode : FileMem::DObject
{
   DNode(): DObject(12)
   {
      val() = 42;
      next() = FileMem::nullPtr;
   }
   uint32_t& val() { get<uint32_t>(0); }
   FileMem::VPtr& next() { get<FileMem::VPtr>(4); }
};



namespace unit
{
   namespace filemem
   {
      constexpr const char* fileName = "testfiles/filememdat";
   }
}



bool unit::filemem::main()
{
   bool ret = false;
   header("FileMem");
   try
   {
      ret = object1()&&
            object1()&&
            object2()&&
            dobject1()&&
            dobject2()&&
            dobject3()&&
            dobject4()&&
            dobject5()&&
            init1()&&
            init2()&&
            sync1()&&
            sync2()&&
            operat1();
   }
   catch (...)
   {
      system("rm -f testfiles/filememdat");
      end();
      throw;
   }
   system("rm -f testfiles/filememdat");
   end();
   return ret;
}



bool unit::filemem::object1()
{
   start();
   Node t;
   bool ret = (*reinterpret_cast<uint32_t*>(&t.bytes[0]))==42;
   return finish(ret,"object1");
}



bool unit::filemem::object2()
{
   start();
   Node t;
   bool ret = t.val()==42;
   return finish(ret,"object2");
}



bool unit::filemem::dobject1()
{
   start();
   DNode t;
   bool ret = t.val()==42;
   return finish(ret,"dobject1");
}



bool unit::filemem::dobject2()
{
   start();
   DNode cpy;
   DNode t(cpy);
   bool ret = t.val()==42;
   return finish(ret,"dobject2");
}



bool unit::filemem::dobject3()
{
   start();
   DNode cpy;
   DNode t(std::move(cpy));
   bool ret = t.val()==42;
   return finish(ret,"dobject3");
}



bool unit::filemem::dobject4()
{
   start();
   DNode cpy;
   DNode t;
   t = cpy;
   bool ret = t.val()==42;
   return finish(ret,"dobject4");
}



bool unit::filemem::dobject5()
{
   start();
   DNode cpy;
   DNode t;
   t = std::move(cpy);
   bool ret = t.val()==42;
   return finish(ret,"dobject5");
}



bool unit::filemem::init1()
{
   start();
   FileMem::Base tf(fileName);
   tf.reserve(4096);
   FileMem::Map<Node> t(tf.allocate(12));
   bool ret = t->val()==42;
   return finish(ret,"init1");
}



bool unit::filemem::init2()
{
   start();
   FileMem::Base tf(fileName);
   FileMem::Map<Node> tmp(tf.allocate(12));
   FileMem::Map<Node> t({FileMem::nullPtr,&tf});
   t = tmp.addr();
   bool ret = t.addr()==tmp.addr();
   return finish(ret,"init2");
}



bool unit::filemem::sync1()
{
   start();
   bool ret = false;
   {
      FileMem::Base tf(fileName);
      tf.clear();
      FileMem::Map<Node> t = tf.allocate(12);
      t->val() = 33;
      t.sync(FileMem::Sync::write);
   }
   {
      FileMem::Base tf(fileName);
      FileMem::Map<Node> t = tf.head();
      t.sync(FileMem::Sync::read);
      ret = t->val()==33;
   }
   return finish(ret,"sync1");
}



bool unit::filemem::sync2()
{
   start();
   bool ret = false;
   {
      FileMem::Base tf(fileName);
      tf.clear();
      FileMem::Map<Node> t = tf.allocate(24);
      t->val() = 33;
      t.sync(FileMem::Sync::write,1);
   }
   {
      FileMem::Base tf(fileName);
      FileMem::Map<Node> t = tf.head();
      t.sync(FileMem::Sync::read,1);
      ret = t->val()==33;
   }
   return finish(ret,"sync2");
}



bool unit::filemem::operat1()
{
   start();
   FileMem::Base tf(fileName);
   FileMem::Map<Node> t = tf.head();
   bool ret = (*t).val()==42;
   return finish(ret,"operat1");
}
