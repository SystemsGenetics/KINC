#include "getopts.h"
#include "unit.h"



namespace unit
{
   namespace getopts
   {
      constexpr auto commStr1 = "one -test=1";
      constexpr auto commStr2 = "-first one --second=1 TWO ---third=1.12 three"
                                " ----fourth=ok -----fifth";
      constexpr auto commStr3 = "no options";
   }
}



bool unit::getopts::main()
{
   bool ret = false;
   header("GetOpts");
   try
   {
      ret = init1()&&
            com_get1()&&
            com_get2()&&
            com_front1()&&
            com_pop1()&&
            com_empty1()&&
            size1()&&
            empty1()&&
            iterate1()&&
            erase1()&&
            key1()&&
            iter_empty1()&&
            operate1()&&
            operate2()&&
            operate3()&&
            operate4()&&
            operate5()&&
            operate6()&&
            operate7()&&
            operate8()&&
            operate9()&&
            operate10()&&
            operate11()&&
            operate12()&&
            operate13()&&
            operate14()&&
            operate15()&&
            operate16()&&
            operate17()&&
            operate18();
   }
   catch (...)
   {
      end();
      throw;
   }
   end();
   return ret;
}



bool unit::getopts::init1()
{
   start();
   GetOpts t(commStr1);
   bool ret = t.com_size()==1;
   return finish(ret,"init1");
}



bool unit::getopts::com_get1()
{
   start();
   GetOpts t(commStr1);
   bool ret = t.com_get({"wrong"})==0;
   return finish(ret,"get_com1");
}



bool unit::getopts::com_get2()
{
   start();
   GetOpts t(commStr1);
   bool ret = t.com_get({"one"})==1;
   return finish(ret,"get_com2");
}



bool unit::getopts::com_front1()
{
   start();
   GetOpts t(commStr1);
   bool ret = t.com_front()==std::string("one");
   return finish(ret,"com_front1");
}



bool unit::getopts::com_pop1()
{
   start();
   GetOpts t(commStr2);
   t.com_pop();
   bool ret = t.com_get({"TWO"})==1;
   return finish(ret,"pop_com1");
}



bool unit::getopts::com_empty1()
{
   start();
   GetOpts t(commStr1);
   t.com_pop();
   bool ret = t.com_empty();
   return finish(ret,"com_empty1");
}



bool unit::getopts::size1()
{
   start();
   GetOpts t(commStr2);
   bool ret = t.size()==5;
   return finish(ret,"size1");
}



bool unit::getopts::empty1()
{
   start();
   GetOpts t(commStr3);
   bool ret = t.empty();
   return finish(ret,"empty1");
}



bool unit::getopts::iterate1()
{
   start();
   GetOpts t(commStr2);
   int count {0};
   for (auto i = t.begin();i!=t.end();++i)
   {
      count++;
   }
   bool ret = count==5;
   return finish(ret,"iterate1");
}



bool unit::getopts::erase1()
{
   start();
   GetOpts t(commStr1);
   t.erase(t.begin());
   bool ret = t.begin()==t.end();
   return finish(ret,"erase1");
}



bool unit::getopts::key1()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   bool ret = i.key()==std::string("first");
   ++i;
   ret = i.key()==std::string("second");
   ++i;
   ret = i.key()==std::string("third");
   ++i;
   ret = i.key()==std::string("fourth");
   ++i;
   ret = i.key()==std::string("fifth");
   return finish(ret,"key1");
}



bool unit::getopts::iter_empty1()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   bool ret = i.empty();
   ++i;
   ret = ret&&!i.empty();
   return finish(ret,"iter_empty1");
}



bool unit::getopts::operate1()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   bool ret = false;
   try
   {
      short n;
      i >> n;
   }
   catch (GetOpts::InvalidType)
   {
      ret = true;
   }
   return finish(ret,"operate1");
}



bool unit::getopts::operate2()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   bool ret = false;
   try
   {
      unsigned short n;
      i >> n;
   }
   catch (GetOpts::InvalidType)
   {
      ret = true;
   }
   return finish(ret,"operate2");
}



bool unit::getopts::operate3()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   bool ret = false;
   try
   {
      int n;
      i >> n;
   }
   catch (GetOpts::InvalidType)
   {
      ret = true;
   }
   return finish(ret,"operate3");
}



bool unit::getopts::operate4()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   bool ret = false;
   try
   {
      unsigned int n;
      i >> n;
   }
   catch (GetOpts::InvalidType)
   {
      ret = true;
   }
   return finish(ret,"operate4");
}



bool unit::getopts::operate5()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   bool ret = false;
   try
   {
      long n;
      i >> n;
   }
   catch (GetOpts::InvalidType)
   {
      ret = true;
   }
   return finish(ret,"operate5");
}



bool unit::getopts::operate6()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   bool ret = false;
   try
   {
      unsigned long n;
      i >> n;
   }
   catch (GetOpts::InvalidType)
   {
      ret = true;
   }
   return finish(ret,"operate6");
}



bool unit::getopts::operate7()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   bool ret = false;
   try
   {
      float n;
      i >> n;
   }
   catch (GetOpts::InvalidType)
   {
      ret = true;
   }
   return finish(ret,"operate7");
}



bool unit::getopts::operate8()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   bool ret = false;
   try
   {
      double n;
      i >> n;
   }
   catch (GetOpts::InvalidType)
   {
      ret = true;
   }
   return finish(ret,"operate8");
}



bool unit::getopts::operate9()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   bool ret = false;
   try
   {
      std::string n;
      i >> n;
   }
   catch (GetOpts::InvalidType)
   {
      ret = true;
   }
   return finish(ret,"operate9");
}



bool unit::getopts::operate10()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   ++i;
   short n;
   i >> n;
   bool ret = n==1;
   return finish(ret,"operate10");
}



bool unit::getopts::operate11()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   ++i;
   unsigned short n;
   i >> n;
   bool ret = n==1;
   return finish(ret,"operate11");
}



bool unit::getopts::operate12()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   ++i;
   int n;
   i >> n;
   bool ret = n==1;
   return finish(ret,"operate12");
}



bool unit::getopts::operate13()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   ++i;
   unsigned int n;
   i >> n;
   bool ret = n==1;
   return finish(ret,"operate13");
}



bool unit::getopts::operate14()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   ++i;
   long n;
   i >> n;
   bool ret = n==1;
   return finish(ret,"operate14");
}



bool unit::getopts::operate15()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   ++i;
   unsigned long n;
   i >> n;
   bool ret = n==1;
   return finish(ret,"operate15");
}



bool unit::getopts::operate16()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   ++i;
   ++i;
   float n;
   i >> n;
   bool ret = n==1.12f;
   return finish(ret,"operate16");
}



bool unit::getopts::operate17()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   ++i;
   ++i;
   double n;
   i >> n;
   bool ret = n==1.12;
   return finish(ret,"operate17");
}



bool unit::getopts::operate18()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   ++i;
   ++i;
   ++i;
   std::string n;
   i >> n;
   bool ret = n==std::string("ok");
   return finish(ret,"operate18");
}
