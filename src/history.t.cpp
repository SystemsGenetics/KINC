#include "history.h"
#include "unit.h"



namespace unit
{
   namespace history
   {
      constexpr auto tmpFile = "testfiles/histfile.tmp";
      constexpr auto tmpFile2 = "testfiles/histfile2.tmp";
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
   std::string rmCmd = "rm -f ";
   rmCmd += tmpFile;
   rmCmd += " ";
   rmCmd += tmpFile2;
   try
   {
      construct_history_file();
      ret = construct()&&
            add_child()&&
            iterate()&&
            iter_childhead();
   }
   catch (...)
   {
      system(rmCmd.c_str());
      unit::end();
      throw;
   }
   system(rmCmd.c_str());
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



bool unit::history::iter_childhead()
{
   start();
   FileMem tf(tmpFile);
   History th(tf,tf.head());
   auto i = th.begin();
   i = i.childHead();
   auto t = i.load();
   bool test = t.timeStamp()==child2Num;
   return finish(test,"iter_childhead1");
}
