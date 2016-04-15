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
      constexpr auto commStr4 = "--no --comms";
      constexpr auto invalidComm = "hello ---test=1=1";
   }
}



bool unit::getopts::main()
{
   bool ret = false;
   header("GetOpts");
   try
   {
      ret = construct()&&
            orig()&&
            com_size()&&
            com_empty()&&
            com_get()&&
            com_front()&&
            com_pop()&&
            com_empty()&&
            size()&&
            empty()&&
            has_opt()&&
            iterate()&&
            erase()&&
            iter_key()&&
            iter_value()&&
            iter_is_key()&&
            iter_val_empty()&&
            iter_operat_equal();
   }
   catch (...)
   {
      end();
      throw;
   }
   end();
   return ret;
}



bool unit::getopts::construct()
{
   bool cont {true};
   {
      start();
      GetOpts t(commStr1);
      bool test = t.com_size()==1;
      cont = cont&&finish(test,"construct1");
   }
   if (cont)
   {
      start();
      bool test {false};
      try
      {
         GetOpts t(invalidComm);
      }
      catch (GetOpts::InvalidSyntax)
      {
         test = true;
      }
      cont = cont&&finish(test,"construct2");
   }
   return cont;
}



bool unit::getopts::orig()
{
   start();
   GetOpts t(commStr1);
   bool test {t.orig()==std::string(commStr1)};
   return finish(test,"orig1");
}



bool unit::getopts::com_size()
{
   start();
   GetOpts t(commStr2);
   bool test {t.com_size()==3};
   return finish(test,"com_size1");
}



bool unit::getopts::com_empty()
{
   start();
   GetOpts t(commStr4);
   bool test {t.com_empty()};
   return finish(test,"com_empty1");
}



bool unit::getopts::com_get()
{
   bool cont {true};
   {
      start();
      GetOpts t(commStr1);
      bool test = t.com_get({"wrong"})==0;
      cont = cont&&finish(test,"get_com1");
   }
   if (cont)
   {
      start();
      GetOpts t(commStr1);
      bool test = t.com_get({"one"})==1;
      cont = cont&&finish(test,"get_com2");
   }
   return cont;
}



bool unit::getopts::com_front()
{
   start();
   GetOpts t(commStr1);
   bool test = t.com_front()==std::string("one");
   return finish(test,"com_front1");
}



bool unit::getopts::com_pop()
{
   start();
   GetOpts t(commStr2);
   t.com_pop();
   bool test = t.com_get({"TWO"})==1;
   return finish(test,"pop_com1");
}



bool unit::getopts::size()
{
   start();
   GetOpts t(commStr2);
   bool test = t.size()==5;
   return finish(test,"size1");
}



bool unit::getopts::empty()
{
   start();
   GetOpts t(commStr3);
   bool test = t.empty();
   return finish(test,"empty1");
}



bool unit::getopts::has_opt()
{
   bool cont {true};
   {
      start();
      GetOpts t(commStr2);
      bool test {t.has_opt("third")};
      cont = cont&&finish(test,"has_opt1");
   }
   if (cont)
   {
      start();
      GetOpts t(commStr2);
      t.has_opt("third",true);
      bool test {t.size()==4};
      cont = cont&&finish(test,"has_opt2");
   }
   return cont;
}



bool unit::getopts::iterate()
{
   start();
   GetOpts t(commStr2);
   int count {0};
   for (auto i = t.begin();i!=t.end();++i)
   {
      count++;
   }
   bool test {count==5};
   return finish(test,"iterate1");
}



bool unit::getopts::erase()
{
   start();
   GetOpts t(commStr1);
   t.erase(t.begin());
   bool test {t.begin()==t.end()};
   return finish(test,"erase1");
}



bool unit::getopts::iter_key()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   bool test = i.key()==std::string("first");
   ++i;
   test = test&&i.key()==std::string("second");
   ++i;
   test = test&&i.key()==std::string("third");
   ++i;
   test = test&&i.key()==std::string("fourth");
   ++i;
   test = test&&i.key()==std::string("fifth");
   return finish(test,"iter_key1");
}



bool unit::getopts::iter_value()
{
   bool cont {true};
   {
      start();
      GetOpts t(commStr2);
      auto i = t.begin();
      bool test {false};
      try
      {
         i.value<int>();
      }
      catch (GetOpts::InvalidType)
      {
         test = true;
      }
      cont = cont&&finish(test,"iter_value1");
   }
   if (cont)
   {
      start();
      GetOpts t(commStr2);
      auto i = t.begin();
      ++i;
      bool test = i.value<int>()==1;
      ++i;
      test = test&&i.value<float>()==1.12f;
      ++i;
      test = test&&i.value<std::string>()==std::string("ok");
      cont = cont&&finish(test,"iter_value2");
   }
   return cont;
}



bool unit::getopts::iter_is_key()
{
   start();
   GetOpts t(commStr1);
   auto i = t.begin();
   bool test {i.is_key("test")};
   return finish(test,"iter_is_key1");
}



bool unit::getopts::iter_val_empty()
{
   start();
   GetOpts t(commStr2);
   auto i = t.begin();
   bool test = i.val_empty();
   ++i;
   test = test&&!i.val_empty();
   return finish(test,"iter_val_empty1");
}



bool unit::getopts::iter_operat_equal()
{
   start();
   GetOpts t(commStr1);
   auto x = t.begin();
   auto y = t.begin();
   bool test {x==y};
   return finish(test,"iter_operat_equal");
}
