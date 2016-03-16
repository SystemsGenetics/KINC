#include "history.h"
#include "unit.h"



namespace unit
{
   namespace history
   {
      constexpr auto tmpFile = "histfile.tmp";
      constexpr auto tmpFile2 = "histfile2.tmp";
      constexpr int headNum = 9999;
      constexpr int childNum = 8888;
      constexpr int child2Num = 7777;
      constexpr int nextNum = 6666;
   }
}



void construct_history_file()
{
   FileMem tf(unit::history::tmpFile);
   HistItem t(tf);
   t.allocate();
   t.timeStamp(unit::history::headNum);
   HistItem ch(tf);
   ch.allocate();
   ch.timeStamp(unit::history::childNum);
   HistItem n(tf);
   n.allocate();
   n.timeStamp(unit::history::nextNum);
   HistItem ch2(tf);
   ch2.allocate();
   ch2.timeStamp(unit::history::child2Num);
   t.childHead(ch.addr());
   ch.next(n.addr());
   ch.childHead(ch2.addr());
   t.sync();
   ch.sync();
   ch2.sync();
   n.sync();
}



bool unit::history::main()
{
   bool ret = false;
   header("History");
   try
   {
      construct_history_file();
      ret = construct()&&
            add_child()&&
            has_child()&&
            iterate()&&
            iter_has_child()&&
            iter_load()&&
            iter_child();
   }
   catch (...)
   {
      system("rm -f *.tmp");
      unit::end();
      throw;
   }
   system("rm -f *.tmp");
   unit::end();
   return ret;
}



bool unit::history::construct()
{
   bool cont = true;
   {
      start();
      FileMem tf(tmpFile);
      History t(tf,tf.head());
      bool test = t.addr()==tf.head();
      cont = cont&&finish(test,"contruct1");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile2);
      History t(tf);
      bool test = t.addr()==tf.head();
      cont = cont&&finish(test,"contruct2");
   }
   return cont;
}



bool unit::history::add_child()
{
   bool cont = true;
   {
      start();
      FileMem tf(tmpFile2);
      History t(tf,tf.head());
      History cp(tf);
      cp.timeStamp(9999);
      t.add_child(cp);
      HistItem th(tf,tf.head());
      th = th.childHead();
      bool test = th.timeStamp()==9999;
      cont = cont&&finish(test,"add_child1");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile2);
      History t(tf,tf.head());
      History cp(tf);
      cp.timeStamp(9999);
      t.add_child(cp);
      HistItem th(tf,tf.head());
      th = th.childHead();
      th = th.next();
      bool test = th.timeStamp()==9999;
      cont = cont&&finish(test,"add_child2");
   }
   return cont;
}



bool unit::history::has_child()
{
   start();
   FileMem tf(tmpFile2);
   tf.clear();
   History t(tf);
   bool test = !t.has_child();
   return finish(test,"has_child1");
}



bool unit::history::iterate()
{
   start();
   FileMem tf(tmpFile);
   History t(tf,tf.head());
   int count = 0;
   for (auto i = t.begin();i!=t.end();++i)
   {
      ++count;
   }
   bool test = count==2;
   return finish(test,"iterate1");
}



bool unit::history::iter_has_child()
{
   start();
   FileMem tf(tmpFile);
   History t(tf,tf.head());
   auto i = t.begin();
   bool test = i.has_child();
   return finish(test,"iter_has_child1");
}



bool unit::history::iter_load()
{
   start();
   FileMem tf(tmpFile);
   History th(tf,tf.head());
   auto i = th.begin();
   auto t = i.load();
   bool test = t.timeStamp()==childNum;
   return finish(test,"iter_load1");
}



bool unit::history::iter_child()
{
   start();
   FileMem tf(tmpFile);
   History th(tf,tf.head());
   auto i = th.begin();
   i = i.child();
   auto t = i.load();
   bool test = t.timeStamp()==child2Num;
   return finish(test,"iter_childhead1");
}
