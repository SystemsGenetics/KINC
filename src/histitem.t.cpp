#include <string>
#include "unit.h"
#include "histitem.h"



namespace unit
{
   namespace histitem
   {
      constexpr auto tmpFile = "testfiles/histfile.tmp";
      constexpr auto tmpFile2 = "testfiles/histfile2.tmp";
      constexpr auto testStr = "1234567890123456";
   }
}



bool unit::histitem::main()
{
   bool ret = false;
   header("HistItem");
   std::string rmCmd = "rm -f ";
   rmCmd += tmpFile;
   rmCmd += " ";
   rmCmd += tmpFile2;
   try
   {
      ret = init1()&&
            allocate1()&&
            init2()&&
            allocate2()&&
            init3()&&
            init4()&&
            timestamp1()&&
            timestamp2()&&
            filename1()&&
            filename2()&&
            filename3()&&
            filename4()&&
            object1()&&
            object2()&&
            object3()&&
            object4()&&
            command1()&&
            command2()&&
            command3()&&
            command4()&&
            next1()&&
            next2()&&
            next3()&&
            next4()&&
            childhead1()&&
            childhead2()&&
            childhead3()&&
            childhead4()&&
            init5();
   }
   catch (...)
   {
      system(rmCmd.c_str());
      end();
      throw;
   }
   system(rmCmd.c_str());
   end();
   return ret;
}



bool unit::histitem::init1()
{
   start();
   FileMem tf(tmpFile);
   HistItem t(tf);
   bool ret = t.addr()==FileMem::nullPtr;
   return finish(ret,"init1");
}



bool unit::histitem::allocate1()
{
   start();
   FileMem tf(tmpFile);
   HistItem t(tf);
   t.allocate();
   bool ret = t.addr()!=FileMem::nullPtr;
   return finish(ret,"allocate1");
}



bool unit::histitem::init2()
{
   start();
   {
      FileMem tf(tmpFile);
      tf.clear();
      HistItem t(tf);
      t.allocate();
      t.sync();
   }
   FileMem tf(tmpFile);
   HistItem t(tf,tf.head());
   bool ret = t.addr()!=FileMem::nullPtr;
   return finish(ret,"init2");
}



bool unit::histitem::allocate2()
{
   start();
   FileMem tf(tmpFile);
   HistItem t(tf,tf.head());
   bool ret = false;
   try
   {
      t.allocate();
   }
   catch (HistItem::IsAllocated)
   {
      ret = true;
   }
   return finish(ret,"allocate2");
}



bool unit::histitem::init3()
{
   start();
   {
      FileMem tf(tmpFile);
      tf.clear();
      HistItem t(tf);
      t.allocate();
      t.timeStamp(9999);
      t.sync();
   }
   FileMem tf(tmpFile);
   HistItem t(tf);
   t = tf.head();
   bool ret = t.timeStamp()==9999;
   return finish(ret,"init3");
}



bool unit::histitem::init4()
{
   start();
   {
      FileMem tf(tmpFile);
      tf.clear();
      HistItem t(tf);
      HistItem ch(tf);
      t.allocate();
      ch.allocate();
      t.timeStamp(9999);
      ch.timeStamp(8888);
      t.next(ch.addr());
      ch.sync();
      ch = FileMem::nullPtr;
      ch.allocate();
      ch.timeStamp(7777);
      t.childHead(ch.addr());
      ch.sync();
      t.sync();
   }
   {
      FileMem tf(tmpFile);
      FileMem tf2(tmpFile2);
      HistItem cp(tf,tf.head());
      HistItem t(tf2);
      t = cp;
   }
   FileMem tf(tmpFile2);
   HistItem t(tf,tf.head());
   bool ret = t.timeStamp()==9999;
   t = t.next();
   ret = ret&&t.timeStamp()==8888;
   t = tf.head();
   t = t.childHead();
   ret = ret&&t.timeStamp()==7777;
   return finish(ret,"init4");
}



bool unit::histitem::timestamp1()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,FileMem::nullPtr);
      t.timeStamp(9999);
   }
   catch (HistItem::IsNullPtr)
   {
      ret = true;
   }
   return finish(ret,"timestamp1");
}



bool unit::histitem::timestamp2()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,FileMem::nullPtr);
      t.timeStamp();
   }
   catch (HistItem::IsNullPtr)
   {
      ret = true;
   }
   return finish(ret,"timestamp2");
}



bool unit::histitem::filename1()
{
   start();
   {
      FileMem tf(tmpFile);
      tf.clear();
      HistItem t(tf);
      t.allocate();
      t.fileName(testStr);
      t.sync();
   }
   FileMem tf(tmpFile);
   HistItem t(tf,tf.head());
   std::string tmp = testStr;
   bool ret = t.fileName()==tmp;
   return finish(ret,"filename1");
}



bool unit::histitem::filename2()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,tf.head());
      t.fileName(testStr);
   }
   catch (HistItem::AlreadySet)
   {
      ret = true;
   }
   return finish(ret,"filename2");
}



bool unit::histitem::filename3()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,FileMem::nullPtr);
      t.fileName(testStr);
   }
   catch (HistItem::IsNullPtr)
   {
      ret = true;
   }
   return finish(ret,"filename3");
}



bool unit::histitem::filename4()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,FileMem::nullPtr);
      t.fileName();
   }
   catch (HistItem::IsNullPtr)
   {
      ret = true;
   }
   return finish(ret,"filename4");
}



bool unit::histitem::object1()
{
   start();
   {
      FileMem tf(tmpFile);
      tf.clear();
      HistItem t(tf);
      t.allocate();
      t.object(testStr);
      t.sync();
   }
   FileMem tf(tmpFile);
   HistItem t(tf,tf.head());
   std::string tmp = testStr;
   bool ret = t.object()==tmp;
   return finish(ret,"object1");
}



bool unit::histitem::object2()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,tf.head());
      t.object(testStr);
   }
   catch (HistItem::AlreadySet)
   {
      ret = true;
   }
   return finish(ret,"object2");
}



bool unit::histitem::object3()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,FileMem::nullPtr);
      t.object(testStr);
   }
   catch (HistItem::IsNullPtr)
   {
      ret = true;
   }
   return finish(ret,"object3");
}



bool unit::histitem::object4()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,FileMem::nullPtr);
      t.object();
   }
   catch (HistItem::IsNullPtr)
   {
      ret = true;
   }
   return finish(ret,"object4");
}



bool unit::histitem::command1()
{
   start();
   {
      FileMem tf(tmpFile);
      tf.clear();
      HistItem t(tf);
      t.allocate();
      t.command(testStr);
      t.sync();
   }
   FileMem tf(tmpFile);
   HistItem t(tf,tf.head());
   std::string tmp = testStr;
   bool ret = t.command()==tmp;
   return finish(ret,"command1");
}



bool unit::histitem::command2()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,tf.head());
      t.command(testStr);
   }
   catch (HistItem::AlreadySet)
   {
      ret = true;
   }
   return finish(ret,"command2");
}



bool unit::histitem::command3()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,FileMem::nullPtr);
      t.command(testStr);
   }
   catch (HistItem::IsNullPtr)
   {
      ret = true;
   }
   return finish(ret,"command3");
}



bool unit::histitem::command4()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,FileMem::nullPtr);
      t.command();
   }
   catch (HistItem::IsNullPtr)
   {
      ret = true;
   }
   return finish(ret,"command4");
}



bool unit::histitem::next1()
{
   start();
   {
      FileMem tf(tmpFile);
      tf.clear();
      HistItem t(tf);
      t.allocate();
      t.next(9999);
      t.sync();
   }
   FileMem tf(tmpFile);
   HistItem t(tf,tf.head());
   bool ret = t.next()==9999;
   return finish(ret,"next1");
}



bool unit::histitem::next2()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,tf.head());
      t.next(9999);
   }
   catch (HistItem::AlreadySet)
   {
      ret = true;
   }
   return finish(ret,"next2");
}



bool unit::histitem::next3()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,FileMem::nullPtr);
      t.next(9999);
   }
   catch (HistItem::IsNullPtr)
   {
      ret = true;
   }
   return finish(ret,"next3");
}



bool unit::histitem::next4()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,FileMem::nullPtr);
      t.next();
   }
   catch (HistItem::IsNullPtr)
   {
      ret = true;
   }
   return finish(ret,"next4");
}



bool unit::histitem::childhead1()
{
   start();
   {
      FileMem tf(tmpFile);
      tf.clear();
      HistItem t(tf);
      t.allocate();
      t.childHead(9999);
      t.sync();
   }
   FileMem tf(tmpFile);
   HistItem t(tf,tf.head());
   bool ret = t.childHead()==9999;
   return finish(ret,"childhead1");
}



bool unit::histitem::childhead2()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,tf.head());
      t.childHead(9999);
   }
   catch (HistItem::AlreadySet)
   {
      ret = true;
   }
   return finish(ret,"childhead2");
}



bool unit::histitem::childhead3()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,FileMem::nullPtr);
      t.childHead(9999);
   }
   catch (HistItem::IsNullPtr)
   {
      ret = true;
   }
   return finish(ret,"childhead3");
}



bool unit::histitem::childhead4()
{
   start();
   bool ret = false;
   try
   {
      FileMem tf(tmpFile);
      HistItem t(tf,FileMem::nullPtr);
      t.childHead();
   }
   catch (HistItem::IsNullPtr)
   {
      ret = true;
   }
   return finish(ret,"childhead4");
}



bool unit::histitem::init5()
{
   start();
   std::string fileName = "filename";
   std::string object = "object";
   std::string command = "command";
   {
      FileMem tf(tmpFile);
      tf.clear();
      HistItem t(tf);
      t.allocate();
      t.timeStamp(9999);
      t.fileName(fileName);
      t.object(object);
      t.command(command);
      t.next(8888);
      t.childHead(7777);
      t.sync();
   }
   FileMem tf(tmpFile);
   HistItem t(tf,tf.head());
   bool ret = t.timeStamp()==9999;
   ret = ret&&t.fileName()==fileName;
   ret = ret&&t.object()==object;
   ret = ret&&t.command()==command;
   ret = ret&&t.next()==8888;
   ret = ret&&t.childHead()==7777;
   return finish(ret,"init5");
}
