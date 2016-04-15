#include <string>
#include "unit.h"
#include "histitem.h"



namespace unit
{
   namespace histitem
   {
      constexpr auto tmpFile = "histfile.tmp";
      constexpr auto tmpFile2 = "histfile2.tmp";
      constexpr auto testStr = "1234567890123456";
   }
}



bool unit::histitem::main()
{
   bool ret = false;
   header("HistItem");
   try
   {
      ret = construct()&&
            move()&&
            allocate()&&
            sync()&&
            timestamp()&&
            filename()&&
            object()&&
            command()&&
            next()&&
            childhead()&&
            operat_set()&&
            copy_from()&&
            mem()&&
            final();
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



bool unit::histitem::construct()
{
   bool cont = true;
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf);
      bool test = t.addr()==FileMem::nullPtr;
      cont = cont&&finish(test,"construct1");
   }
   if (cont)
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
      bool test = t.addr()!=FileMem::nullPtr;
      cont = cont&&finish(test,"construct2");
   }
   return cont;
}



bool unit::histitem::move()
{
   bool cont = true;
   {
      start();
      FileMem tf(tmpFile);
      HistItem tmp(tf,tf.head());
      HistItem t(std::move(tmp));
      bool test = t.addr()!=FileMem::nullPtr&&tmp.addr()==FileMem::nullPtr;
      cont = cont&&finish(test,"move1");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem tmp(tf,tf.head());
      HistItem t(tf);
      t = std::move(tmp);
      bool test = t.addr()!=FileMem::nullPtr&&tmp.addr()==FileMem::nullPtr;
      cont = cont&&finish(test,"move2");
   }
   return cont;
}



bool unit::histitem::allocate()
{
   bool cont = true;
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf);
      t.allocate();
      bool test = t.addr()!=FileMem::nullPtr;
      cont = cont&&finish(test,"allocate1");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf);
      t.allocate();
      bool test = false;
      try
      {
         t.allocate();
      }
      catch (HistItem::IsAllocated)
      {
         test = true;
      }
      cont = cont&&finish(test,"allocate2");
   }
   return cont;
}



bool unit::histitem::sync()
{
   bool cont = true;
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
      HistItem t(tf,tf.head());
      bool test = t.timeStamp()==9999;
      cont = cont&&finish(test,"sync1");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf);
      bool test = false;
      try
      {
         t.sync();
      }
      catch (HistItem::IsNullPtr)
      {
         test = true;
      }
      cont = cont&&finish(test,"sync2");
   }
   return cont;
}



bool unit::histitem::timestamp()
{
   start();
   FileMem tf(tmpFile);
   HistItem t(tf);
   bool test = false;
   try
   {
      t.timeStamp(9999);
   }
   catch (HistItem::IsNullPtr)
   {
      test = true;
   }
   return finish(test,"timestamp1");
}



bool unit::histitem::filename()
{
   bool cont = true;
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
      bool test = strcmp(t.fileName().c_str(),testStr)==0;
      cont = cont&&finish(test,"filename1");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf,tf.head());
      bool test = false;
      try
      {
         t.fileName(testStr);
      }
      catch (HistItem::AlreadySet)
      {
         test = true;
      }
      cont = cont&&finish(test,"filename2");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf);
      bool test = false;
      try
      {
         t.fileName(testStr);
      }
      catch (HistItem::IsNullPtr)
      {
         test = true;
      }
      cont = cont&&finish(test,"filename3");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf);
      bool test = false;
      try
      {
         t.fileName();
      }
      catch (HistItem::IsNullPtr)
      {
         test = true;
      }
      cont = cont&&finish(test,"filename4");
   }
   return cont;
}



bool unit::histitem::object()
{
   bool cont = true;
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
      bool test = strcmp(t.object().c_str(),testStr)==0;
      cont = cont&&finish(test,"object1");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf,tf.head());
      bool test = false;
      try
      {
         t.object(testStr);
      }
      catch (HistItem::AlreadySet)
      {
         test = true;
      }
      cont = cont&&finish(test,"object2");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf);
      bool test = false;
      try
      {
         t.object(testStr);
      }
      catch (HistItem::IsNullPtr)
      {
         test = true;
      }
      cont = cont&&finish(test,"object3");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf);
      bool test = false;
      try
      {
         t.object();
      }
      catch (HistItem::IsNullPtr)
      {
         test = true;
      }
      cont = cont&&finish(test,"object4");
   }
   return cont;
}



bool unit::histitem::command()
{
   bool cont = true;
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
      bool test = strcmp(t.command().c_str(),testStr)==0;
      cont = cont&&finish(test,"command1");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf,tf.head());
      bool test = false;
      try
      {
         t.command(testStr);
      }
      catch (HistItem::AlreadySet)
      {
         test = true;
      }
      cont = cont&&finish(test,"command2");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf);
      bool test = false;
      try
      {
         t.command(testStr);
      }
      catch (HistItem::IsNullPtr)
      {
         test = true;
      }
      cont = cont&&finish(test,"command3");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf);
      bool test = false;
      try
      {
         t.command();
      }
      catch (HistItem::IsNullPtr)
      {
         test = true;
      }
      cont = cont&&finish(test,"command4");
   }
   return cont;
}



bool unit::histitem::next()
{
   bool cont = true;
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
      bool test = t.next()==9999;
      cont = cont&&finish(test,"next1");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf,tf.head());
      bool test = false;
      try
      {
         t.next(9999);
      }
      catch (HistItem::AlreadySet)
      {
         test = true;
      }
      cont = cont&&finish(test,"next2");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf);
      bool test = false;
      try
      {
         t.next(9999);
      }
      catch (HistItem::IsNullPtr)
      {
         test = true;
      }
      cont = cont&&finish(test,"next3");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf);
      bool test = false;
      try
      {
         t.next();
      }
      catch (HistItem::IsNullPtr)
      {
         test = true;
      }
      cont = cont&&finish(test,"next4");
   }
   return cont;
}



bool unit::histitem::childhead()
{
   bool cont = true;
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
      bool test = t.childHead()==9999;
      cont = cont&&finish(test,"childhead1");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf,tf.head());
      bool test = false;
      try
      {
         t.childHead(9999);
      }
      catch (HistItem::AlreadySet)
      {
         test = true;
      }
      cont = cont&&finish(test,"childhead2");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf);
      bool test = false;
      try
      {
         t.childHead(9999);
      }
      catch (HistItem::IsNullPtr)
      {
         test = true;
      }
      cont = cont&&finish(test,"childhead3");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf);
      bool test = false;
      try
      {
         t.childHead();
      }
      catch (HistItem::IsNullPtr)
      {
         test = true;
      }
      cont = cont&&finish(test,"childhead4");
   }
   return cont;
}



bool unit::histitem::operat_set()
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
   HistItem t(tf);
   t = tf.head();
   bool test = t.childHead()==9999;
   return finish(test,"operat_set1");
}



bool unit::histitem::copy_from()
{
   bool cont {true};
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
         t.copy_from(cp);
      }
      FileMem tf(tmpFile2);
      HistItem t(tf,tf.head());
      bool test = t.timeStamp()==9999;
      t = t.next();
      test = test&&t.timeStamp()==8888;
      t = tf.head();
      t = t.childHead();
      test = test&&t.timeStamp()==7777;
      cont = cont&&finish(test,"copy_from1");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf,tf.head());
      HistItem cp(tf,t.next());
      bool test = false;
      try
      {
         t.copy_from(cp);
      }
      catch (HistItem::IsAllocated)
      {
         test = true;
      }
      cont = cont&&finish(test,"copy_from2");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      HistItem t(tf);
      HistItem cp(tf);
      bool test = false;
      try
      {
         t.copy_from(cp);
      }
      catch (HistItem::IsNullPtr)
      {
         test = true;
      }
      cont = cont&&finish(test,"copy_from3");
   }
   return cont;
}



bool unit::histitem::mem()
{
   start();
   FileMem tf(tmpFile);
   HistItem t(tf);
   bool test = t.mem()==&tf;
   return finish(test,"mem1");
}



bool unit::histitem::final()
{
   bool cont = true;
   std::string fileName = "filename";
   std::string object = "object";
   std::string command = "command";
   {
      start();
      FileMem tf(tmpFile);
      tf.clear();
      HistItem tmp(tf);
      tmp.allocate();
      tmp.timeStamp(9999);
      tmp.fileName(fileName);
      tmp.object(object);
      tmp.command(command);
      tmp.next(8888);
      tmp.childHead(7777);
      tmp.sync();
      HistItem t(std::move(tmp));
      bool test = t.timeStamp()==9999;
      test = test&&t.fileName()==fileName;
      test = test&&t.object()==object;
      test = test&&t.command()==command;
      test = test&&t.next()==8888;
      test = test&&t.childHead()==7777;
      test = test&&t.addr()!=FileMem::nullPtr&&tmp.addr()==FileMem::nullPtr;
      cont = cont&&finish(test,"final1");
   }
   if (cont)
   {
      start();
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
      bool test = t.timeStamp()==9999;
      test = test&&t.fileName()==fileName;
      test = test&&t.object()==object;
      test = test&&t.command()==command;
      test = test&&t.next()==8888;
      test = test&&t.childHead()==7777;
      cont = cont&&finish(test,"final2");
   }
   if (cont)
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
         t.fileName(fileName);
         t.object(object);
         t.command(command);
         ch.timeStamp(8888);
         ch.fileName(fileName);
         ch.object(object);
         ch.command(command);
         t.next(ch.addr());
         ch.sync();
         ch = FileMem::nullPtr;
         ch.allocate();
         ch.timeStamp(7777);
         ch.fileName(fileName);
         ch.object(object);
         ch.command(command);
         t.childHead(ch.addr());
         ch.sync();
         t.sync();
      }
      {
         FileMem tf(tmpFile);
         FileMem tf2(tmpFile2);
         tf2.clear();
         HistItem cp(tf,tf.head());
         HistItem t(tf2);
         t.copy_from(cp);
      }
      FileMem tf(tmpFile2);
      HistItem t(tf,tf.head());
      bool test = t.timeStamp()==9999;
      test = test&&t.fileName()==fileName;
      test = test&&t.object()==object;
      test = test&&t.command()==command;
      t = t.next();
      test = test&&t.timeStamp()==8888;
      test = test&&t.fileName()==fileName;
      test = test&&t.object()==object;
      test = test&&t.command()==command;
      t = tf.head();
      t = t.childHead();
      test = test&&t.timeStamp()==7777;
      test = test&&t.fileName()==fileName;
      test = test&&t.object()==object;
      test = test&&t.command()==command;
      cont = cont&&finish(test,"final3");
   }
   return cont;
}
