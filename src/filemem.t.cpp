#include <cstdlib>
#include <iostream>
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



bool unit::filemem::main()
{
   bool ret = false;
   std::cout << "FileMem";
   try
   {
      ret = init1()&&
            init2()&&
            init3()&&
            init4()&&
            init5()&&
            opmove1()&&
            opref1()&&
            opref2()&&
            opcallref1()&&
            opcallref2()&&
            node1()&&
            node2()&&
            node3()&&
            node4()&&
            node5()&&
            node6()&&
            node7()&&
            node8()&&
            node9()&&
            node10();
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



bool unit::filemem::init1()
{
   std::cout << "." << std::flush;
   bool ret = false;
   LinuxFile tf("testfiles/filememdat");
   tf.reserve(4096);
   FileMem::Ptr<Node> t(tf);
   ret = t->val()==42;
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::init1 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::init2()
{
   std::cout << "." << std::flush;
   bool ret = false;
   LinuxFile tf("testfiles/filememdat");
   tf.reserve(4096);
   Node tmp;
   tmp.val() = 43;
   FileMem::Ptr<Node> t(tf,tmp);
   ret = t->val()==43;
   t.save();
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::init2 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::init3()
{
   std::cout << "." << std::flush;
   bool ret = false;
   LinuxFile tf("testfiles/filememdat");
   Node tmp;
   tmp.val() = 43;
   FileMem::Ptr<Node> t(tf,std::move(tmp));
   ret = t->val()==43;
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::init3 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::init4()
{
   std::cout << "." << std::flush;
   bool ret = false;
   LinuxFile tf("testfiles/filememdat");
   FileMem::Ptr<Node> t(tf,tf.head()+12);
   ret = t->val()==43;
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::init4 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::init5()
{
   std::cout << "." << std::flush;
   bool ret = false;
   LinuxFile tf("testfiles/filememdat");
   FileMem::Ptr<Node> tmp(tf,tf.head()+12);
   FileMem::Ptr<Node> t(std::move(tmp));
   ret = t->val()==43;
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::init4 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::opmove1()
{
   std::cout << "." << std::flush;
   bool ret = false;
   LinuxFile tf("testfiles/filememdat");
   FileMem::Ptr<Node> tmp(tf,tf.head()+12);
   FileMem::Ptr<Node> t(tf);
   t = std::move(tmp);
   ret = t->val()==43;
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::init4 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::opref1()
{
   std::cout << "." << std::flush;
   bool ret = false;
   LinuxFile tf("testfiles/filememdat");
   FileMem::Ptr<Node> t(tf,tf.head()+12);
   ret = (*t).val()==43;
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::init4 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::opref2()
{
   std::cout << "." << std::flush;
   bool ret = false;
   LinuxFile tf("testfiles/filememdat");
   FileMem::Ptr<Node> t(tf,tf.head());
   (*t).val() = 55;
   ret = (*t).val()==55;
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::init4 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::opcallref1()
{
   std::cout << "." << std::flush;
   bool ret = false;
   LinuxFile tf("testfiles/filememdat");
   FileMem::Ptr<Node> t(tf,tf.head()+12);
   ret = t->val()==43;
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::init4 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::opcallref2()
{
   std::cout << "." << std::flush;
   bool ret = false;
   LinuxFile tf("testfiles/filememdat");
   FileMem::Ptr<Node> t(tf,tf.head());
   t->val() = 55;
   ret = t->val()==55;
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::init4 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::node1()
{
   std::cout << "." << std::flush;
   bool ret = false;
   {
      LinuxFile tf("testfiles/filememdat");
      tf.clear();
      FileMem::Ptr<DNode> t(tf);
      t->val() = 0;
      t->next() = FileMem::nullPtr;
      t.save();
   }
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<DNode> t(tf,tf.head());
      ret = t->val()==0&&t->next()==FileMem::nullPtr;
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::node1 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::node2()
{
   std::cout << "." << std::flush;
   bool ret = false;
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<DNode> link(tf,tf.head());
      while (link->next()!=FileMem::nullPtr)
      {
         link = link->next();
      }
      FileMem::Ptr<DNode> t(tf);
      link->next() = t.addr();
      link.save();
      t->val() = 1;
      t->next() = FileMem::nullPtr;
      t.save();
   }
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<DNode> t(tf,tf.head());
      ret = true;
      for (int i=0;i<2;i++)
      {
         ret = ret&&t->val()==i;
         t = t->next();
      }
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::node2 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::node3()
{
   std::cout << "." << std::flush;
   bool ret = false;
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<DNode> link(tf,tf.head());
      while (link->next()!=FileMem::nullPtr)
      {
         link = link->next();
      }
      FileMem::Ptr<DNode> t(tf);
      link->next() = t.addr();
      link.save();
      t->val() = 2;
      t->next() = FileMem::nullPtr;
      t.save();
   }
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<DNode> t(tf,tf.head());
      ret = true;
      for (int i=0;i<3;i++)
      {
         ret = ret&&t->val()==i;
         t = t->next();
      }
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::node3 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::node4()
{
   std::cout << "." << std::flush;
   bool ret = false;
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<DNode> link(tf,tf.head());
      while (link->next()!=FileMem::nullPtr)
      {
         link = link->next();
      }
      FileMem::Ptr<DNode> t(tf);
      link->next() = t.addr();
      link.save();
      t->val() = 3;
      t->next() = FileMem::nullPtr;
      t.save();
   }
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<DNode> t(tf,tf.head());
      ret = true;
      for (int i=0;i<4;i++)
      {
         ret = ret&&t->val()==i;
         t = t->next();
      }
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::node4 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::node5()
{
   std::cout << "." << std::flush;
   bool ret = false;
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<DNode> link(tf,tf.head());
      while (link->next()!=FileMem::nullPtr)
      {
         link = link->next();
      }
      FileMem::Ptr<DNode> t(tf);
      link->next() = t.addr();
      link.save();
      t->val() = 4;
      t->next() = FileMem::nullPtr;
      t.save();
   }
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<DNode> t(tf,tf.head());
      ret = true;
      for (int i=0;i<5;i++)
      {
         ret = ret&&t->val()==i;
         t = t->next();
      }
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::node5 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::node6()
{
   std::cout << "." << std::flush;
   bool ret = false;
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<Node> link(tf,tf.head());
      while (link->next()!=FileMem::nullPtr)
      {
         link = link->next();
      }
      FileMem::Ptr<Node> t(tf);
      link->next() = t.addr();
      link.save();
      t->val() = 5;
      t->next() = FileMem::nullPtr;
      t.save();
   }
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<Node> t(tf,tf.head());
      ret = true;
      for (int i=0;i<6;i++)
      {
         ret = ret&&t->val()==i;
         t = t->next();
      }
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::node6 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::node7()
{
   std::cout << "." << std::flush;
   bool ret = false;
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<Node> link(tf,tf.head());
      while (link->next()!=FileMem::nullPtr)
      {
         link = link->next();
      }
      FileMem::Ptr<Node> t(tf);
      link->next() = t.addr();
      link.save();
      t->val() = 6;
      t->next() = FileMem::nullPtr;
      t.save();
   }
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<Node> t(tf,tf.head());
      ret = true;
      for (int i=0;i<7;i++)
      {
         ret = ret&&t->val()==i;
         t = t->next();
      }
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::node7 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::node8()
{
   std::cout << "." << std::flush;
   bool ret = false;
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<Node> link(tf,tf.head());
      while (link->next()!=FileMem::nullPtr)
      {
         link = link->next();
      }
      FileMem::Ptr<Node> t(tf);
      link->next() = t.addr();
      link.save();
      t->val() = 7;
      t->next() = FileMem::nullPtr;
      t.save();
   }
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<Node> t(tf,tf.head());
      ret = true;
      for (int i=0;i<8;i++)
      {
         ret = ret&&t->val()==i;
         t = t->next();
      }
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::node8 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::node9()
{
   std::cout << "." << std::flush;
   bool ret = false;
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<Node> link(tf,tf.head());
      while (link->next()!=FileMem::nullPtr)
      {
         link = link->next();
      }
      FileMem::Ptr<Node> t(tf);
      link->next() = t.addr();
      link.save();
      t->val() = 8;
      t->next() = FileMem::nullPtr;
      t.save();
   }
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<Node> t(tf,tf.head());
      ret = true;
      for (int i=0;i<9;i++)
      {
         ret = ret&&t->val()==i;
         t = t->next();
      }
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::node9 FAILED." << std::endl;
   }
   return ret;
}



bool unit::filemem::node10()
{
   std::cout << "." << std::flush;
   bool ret = false;
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<Node> link(tf,tf.head());
      while (link->next()!=FileMem::nullPtr)
      {
         link = link->next();
      }
      FileMem::Ptr<Node> t(tf);
      link->next() = t.addr();
      link.save();
      t->val() = 9;
      t->next() = FileMem::nullPtr;
      t.save();
   }
   {
      LinuxFile tf("testfiles/filememdat");
      FileMem::Ptr<Node> t(tf,tf.head());
      ret = true;
      for (int i=0;i<10;i++)
      {
         ret = ret&&t->val()==i;
         t = t->next();
      }
   }
   if (!ret)
   {
      std::cout << std::endl << "unit::filemem::node10 FAILED." << std::endl;
   }
   return ret;
}
