#include <cstdlib>
#include <utility>
#include "linuxfile.h"
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
            init3()&&
            init4()&&
            init5()&&
            init6()&&
            addr1()&&
            raw1()&&
            raw2()&&
            save1()&&
            operat1()&&
            operat2()&&
            operat3()&&
            operat4()&&
            operat5()&&
            operat6()&&
            operat7()&&
            operat8()&&
            operat9();
   }
   catch (...)
   {
      system("rm -f testfiles/filememdat");
      std::cout << std::endl;
      throw;
   }
   system("rm -f testfiles/filememdat");
   std::cout << std::endl;
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
   LinuxFile tf(fileName);
   tf.reserve(4096);
   FileMem::Ptr<Node> tmp(tf);
   FileMem::Ptr<Node> t(std::move(tmp));
   bool ret = t->val()==42;
   t.save();
   return finish(ret,"init1");
}



bool unit::filemem::init2()
{
   start();
   LinuxFile tf(fileName);
   FileMem::Ptr<Node> tmp(tf);
   FileMem::Ptr<Node> t(tf);
   t = std::move(tmp);
   bool ret = t->val()==42;
   return finish(ret,"init2");
}



bool unit::filemem::init3()
{
   start();
   LinuxFile tf(fileName);
   FileMem::Ptr<Node> t(tf);
   t = tf.head();
   bool ret = t->val()==42;
   return finish(ret,"init3");
}



bool unit::filemem::init4()
{
   start();
   LinuxFile tf(fileName);
   FileMem::Ptr<Node> t(tf);
   bool ret = t->val()==42;
   t.save();
   return finish(ret,"init4");
}



bool unit::filemem::init5()
{
   start();
   LinuxFile tf(fileName);
   FileMem::Ptr<Node> t(tf,tf.head());
   bool ret = t->val()==42;
   return finish(ret,"init5");
}



bool unit::filemem::init6()
{
   start();
   LinuxFile tf(fileName);
   bool ret = false;
   try
   {
      FileMem::Ptr<Node> t(tf,tf.head(),0);
   }
   catch (FileMem::FileSegFault)
   {
      ret = true;
   }
   return finish(ret,"init6");
}



bool unit::filemem::addr1()
{
   start();
   LinuxFile tf(fileName);
   FileMem::Ptr<Node> t(tf,tf.head());
   bool ret = t.addr()==tf.head();
   return finish(ret,"addr1");
}



bool unit::filemem::raw1()
{
   start();
   LinuxFile tf(fileName);
   FileMem::Ptr<Node> t(tf,tf.head(),5);
   t.raw(3);
   bool ret = t->val()==42;
   return finish(ret,"raw1");
}



bool unit::filemem::raw2()
{
   start();
   LinuxFile tf(fileName);
   bool ret = false;
   try
   {
      FileMem::Ptr<Node> t(tf,tf.head());
      t.raw(1);
   }
   catch (FileMem::FileSegFault)
   {
      ret = true;
   }
   return finish(ret,"raw2");
}



bool unit::filemem::save1()
{
   start();
   LinuxFile tf(fileName);
   FileMem::VPtr ptr = FileMem::nullPtr;
   {
      FileMem::Ptr<Node> t(tf);
      ptr = t.addr();
      t.save();
   }
   FileMem::Ptr<Node> t(tf,ptr);
   bool ret = t->val()==42;
   return finish(ret,"save1");
}



bool unit::filemem::operat1()
{
   start();
   LinuxFile tf(fileName);
   FileMem::Ptr<Node> t(tf);
   bool ret = (*t).val()==42;
   return finish(ret,"operat1");
}



bool unit::filemem::operat2()
{
   start();
   LinuxFile tf(fileName);
   FileMem::Ptr<Node> t(tf);
   bool ret = t->val()==42;
   return finish(ret,"operat2");
}



bool unit::filemem::operat3()
{
   start();
   LinuxFile tf(fileName);
   FileMem::Ptr<Node> t(tf);
   bool ret = t[0].val()==42;
   return finish(ret,"operat3");
}



bool unit::filemem::operat4()
{
   start();
   LinuxFile tf(fileName);
   FileMem::VPtr ptr = FileMem::nullPtr;
   {
      FileMem::Ptr<Node> t(tf,FileMem::nullPtr,3);
      t[2].val() = 33;
      t.save();
      ptr = t.addr();
   }
   FileMem::Ptr<Node> t(tf,ptr,3);
   bool ret = t[2].val()==33;
   return finish(ret,"operat4");
}



bool unit::filemem::operat5()
{
   start();
   LinuxFile tf(fileName);
   bool ret = false;
   try
   {
      FileMem::Ptr<Node> t(tf,tf.head());
      t[1].next();
   }
   catch (FileMem::FileSegFault)
   {
      ret = true;
   }
   return finish(ret,"operat5");
}



bool unit::filemem::operat6()
{
   start();
   LinuxFile tf(fileName);
   FileMem::VPtr ptr = FileMem::nullPtr;
   {
      FileMem::Ptr<Node> t(tf);
      t->val() = 33;
      t.save();
      ptr = t.addr();
   }
   FileMem::Ptr<Node> t(tf,ptr);
   ++t;
   bool ret = t->val()==33;
   return finish(ret,"operat6");
}



bool unit::filemem::operat7()
{
   start();
   LinuxFile tf(fileName);
   FileMem::VPtr ptr = FileMem::nullPtr;
   {
      FileMem::Ptr<Node> t(tf,FileMem::nullPtr,2);
      t[1].val() = 33;
      t.save();
      ptr = t.addr();
   }
   FileMem::Ptr<Node> t(tf,ptr,2);
   ++t;
   bool ret = t->val()==33;
   return finish(ret,"operat7");
}



bool unit::filemem::operat8()
{
   start();
   LinuxFile tf(fileName);
   FileMem::VPtr ptr = FileMem::nullPtr;
   {
      FileMem::Ptr<Node> t(tf);
      t->val() = 33;
      t.save();
      ptr = t.addr();
   }
   FileMem::Ptr<Node> t(tf,ptr);
   --t;
   bool ret = t->val()==33;
   return finish(ret,"operat8");
}



bool unit::filemem::operat9()
{
   start();
   LinuxFile tf(fileName);
   FileMem::VPtr ptr = FileMem::nullPtr;
   {
      FileMem::Ptr<Node> t(tf);
      t->val() = 33;
      t.save();
      ptr = t.addr();
   }
   FileMem::Ptr<Node> t(tf,ptr,2);
   t.raw(1);
   --t;
   bool ret = t->val()==33;
   return finish(ret,"operat9");
}
