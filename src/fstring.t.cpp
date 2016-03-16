#include "fstring.h"
#include "unit.h"



namespace unit
{
   namespace fstring
   {
      constexpr const char* testStr = "hello world! A really long string.";
      struct TestString : FileMem::Static<38>
      {
         TestString()
         {
            stripe() = FStringData::strip;
            sSize() = 35;
            memcpy(c_str(),testStr,35);
         }
         uint8_t& stripe() { get<uint8_t>(0); }
         uint16_t& sSize() { get<uint16_t>(1); }
         char* c_str() { &get<char>(3); }
      };
      constexpr auto tmpFile = "strfile.tmp";
   }
}



bool unit::fstring::main()
{
   bool ret = false;
   header("FString");
   try
   {
      ret = construct()&&
            move()&&
            addr()&&
            operat_fp()&&
            operat_set()&&
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



bool unit::fstring::construct()
{
   std::string empty;
   std::string hello(testStr);
   bool cont = true;
   {
      start();
      FileMem tf(tmpFile);
      FString t(&tf);
      bool test = *t==empty;
      cont = cont&&finish(test,"construct1");
   }
   if (cont)
   {
      start();
      {
         FileMem tf(tmpFile);
         tf.clear();
         TestString tmp;
         tf.allot(tmp);
         tf.sync(tmp,FileSync::write);
      }
      FileMem tf(tmpFile);
      FString t(&tf,tf.head());
      bool test = *t==hello;
      cont = cont&&finish(test,"construct2");
   }
   if (cont)
   {
      start();
      {
         FileMem tf(tmpFile);
         tf.clear();
         TestString tmp;
         tmp.stripe() = FStringData::strip+1;
         tf.allot(tmp);
         tf.sync(tmp,FileSync::write);
      }
      bool test = false;
      try
      {
         FileMem tf(tmpFile);
         FString t(&tf,tf.head());
      }
      catch (FString::InvalidPtr)
      {
         test = true;
      }
      cont = cont&&finish(test,"construct3");
   }
   if (cont)
   {
      start();
      bool test = false;
      try
      {
         FString t(nullptr);
      }
      catch(FString::InvalidPtr)
      {
         test = true;
      }
      cont = cont&&finish(test,"construct4");
   }
   return cont;
}



bool unit::fstring::move()
{
   std::string hello(testStr);
   bool cont = true;
   {
      start();
      {
         FileMem tf(tmpFile);
         tf.clear();
         TestString tmp;
         tf.allot(tmp);
         tf.sync(tmp,FileSync::write);
      }
      FileMem tf(tmpFile);
      FString tmp(&tf,tf.head());
      FString t(std::move(tmp));
      bool test = *t==hello;
      cont = cont&&finish(test,"move1");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      FString tmp(&tf,tf.head());
      FString t(&tf);
      t = std::move(tmp);
      bool test = *t==hello;
      cont = cont&&finish(test,"move2");
   }
   return cont;
}



bool unit::fstring::addr()
{
   start();
   {
      FileMem tf(tmpFile);
      tf.clear();
      TestString tmp;
      tf.allot(tmp);
      tf.sync(tmp,FileSync::write);
   }
   FileMem tf(tmpFile);
   FString t(&tf);
   t.addr(tf.head());
   bool test = t.addr()==tf.head();
   return finish(test,"addr1");
}



bool unit::fstring::operat_fp()
{
   start();
   FileMem tf(tmpFile);
   FString t(&tf);
   bool test = t->empty();
   return finish(test,"operat_fp1");
}



bool unit::fstring::operat_set()
{
   std::string hello(testStr);
   bool cont = true;
   {
      start();
      FileMem tf(tmpFile);
      tf.clear();
      FString t(&tf);
      t = testStr;
      bool test = *t==hello;
      cont = cont&&finish(test,"operat_set1");
   }
   if (cont)
   {
      start();
      FileMem tf(tmpFile);
      FString t(&tf,tf.head());
      bool test = false;
      try
      {
         t = testStr;
      }
      catch (FString::AlreadySet)
      {
         test = true;
      }
      cont = cont&&finish(test,"operat_set2");
   }
   return cont;
}



bool unit::fstring::final()
{
   std::string hello(testStr);
   start();
   FileMem::Ptr adr1;
   FileMem::Ptr adr2;
   FileMem::Ptr adr3;
   {
      FileMem tf(tmpFile);
      FString t(&tf);
      t = testStr;
      adr1 = t.addr();
      t.addr(FileMem::nullPtr);
      t = testStr;
      adr2 = t.addr();
      t.addr(FileMem::nullPtr);
      t = testStr;
      adr3 = t.addr();
   }
   FileMem tf(tmpFile);
   FString t(&tf,adr1);
   bool test = *t==hello;
   t.addr(adr2);
   test = test&&*t==hello;
   t.addr(adr3);
   test = test&&*t==hello;
   return finish(test,"final1");
}
