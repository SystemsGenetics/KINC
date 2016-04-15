#include "datamap.h"
#include "unit.h"
#include "plugins/plugins.h"
#include "linuxterm.h"



class FakeData : public DataPlugin
{
public:
   using string = std::string;
   FakeData(const string& type, const string& name):
      DataPlugin(type,name)
   {}
   void load(GetOpts&,Terminal&) override final
   {
      if (touched)
      {
         throw int(0);
      }
      touched = true;
   }
   void dump(GetOpts&,Terminal&) override final
   {
      if (touched)
      {
         throw int(1);
      }
      touched = true;
   }
   void query(GetOpts&,Terminal&) override final
   {
      if (touched)
      {
         throw int(2);
      }
      touched = true;
   }
   bool empty() override final { return true; }
   bool touched {false};
};



DataPlugin* KINCPlugins::new_data(const std::string& type,
                                  const std::string& name)
{
   DataPlugin* ret = nullptr;
   if (strcmp(type.c_str(),"FakeData")==0) ret = new FakeData(type,name);
   return ret;
}



bool unit::datamap::main()
{
   bool ret = false;
   header("DataMap");
   LinuxTerm::stty_raw();
   try
   {
      ret = construct()&&
            open()&&
            close()&&
            unselect()&&
            load()&&
            dump()&&
            query()&&
            iterate()&&
            selected()&&
            iter_file()&&
            iter_type();
   }
   catch (...)
   {
      system("rm -fr *.tmp");
      LinuxTerm::stty_cooked();
      end();
      throw;
   }
   system("rm -fr *.tmp");
   LinuxTerm::stty_cooked();
   end();
   return ret;
}



bool unit::datamap::construct()
{
   start();
   bool test = false;
   try
   {
      DataMap t;
      DataMap tt;
   }
   catch (DataMap::InvalidUse)
   {
      test = true;
   }
   return finish(test,"construct1");
}



bool unit::datamap::open()
{
   DataMap t;
   bool cont = true;
   {
      start();
      DataPlugin* n = t.open("test.tmp","FakeData",true);
      bool test = n==t.current();
      cont = cont&&finish(test,"open1");
   }
   if (cont)
   {
      start();
      bool test = false;
      try
      {
         t.open("test2.tmp","UnknownData");
      }
      catch (DataMap::InvalidType)
      {
         test = true;
      }
      cont = cont&&finish(test,"open2");
   }
   if (cont)
   {
      start();
      bool test = false;
      try
      {
         t.open("test.tmp","FakeData");
      }
      catch (DataMap::AlreadyExists)
      {
         test = true;
      }
      cont = cont&&finish(test,"open3");
   }
   return cont;
}



bool unit::datamap::close()
{
   bool cont {true};
   DataMap t;
   {
      start();
      t.open("test.tmp","FakeData");
      t.close("test.tmp");
      bool test = t.find("test.tmp")==nullptr;
      cont = cont&&finish(test,"close1");
   }
   if (cont)
   {
      start();
      t.open("test.tmp","FakeData",true);
      bool test = t.close("test.tmp");
      cont = cont&&finish(test,"close2");
   }
   return cont;
}



bool unit::datamap::unselect()
{
   DataMap t;
   bool cont {true};
   {
      start();
      t.open("test.tmp","FakeData",true);
      bool test = t.unselect();
      cont = cont&&finish(test,"unselect1");
   }
   if (cont)
   {
      start();
      bool test = !t.unselect();
      cont = cont&&finish(test,"unselect2");
   }
}



bool unit::datamap::load()
{
   LinuxTerm term;
   GetOpts got("one");
   DataMap t;
   t.open("test.tmp","FakeData");
   bool cont = true;
   {
      start();
      bool test = false;
      try
      {
         t.load(got,term);
      }
      catch (DataMap::NoSelect)
      {
         test = true;
      }
      cont = cont&&finish(test,"load1");
   }
   if (cont)
   {
      start();
      t.select("test.tmp");
      t.load(got,term);
      FakeData* tmp = dynamic_cast<FakeData*>(t.find("test.tmp"));
      bool test = tmp->touched;
      cont = cont&&finish(test,"load2");
   }
   if (cont)
   {
      start();
      bool test = false;
      try
      {
         t.load(got,term);
      }
      catch (int x)
      {
         test = x==0;
      }
      cont = cont&&finish(test,"load3");
   }
   return cont;
}



bool unit::datamap::dump()
{
   LinuxTerm term;
   GetOpts got("one");
   DataMap t;
   t.open("test.tmp","FakeData");
   bool cont = true;
   {
      start();
      bool test = false;
      try
      {
         t.dump(got,term);
      }
      catch (DataMap::NoSelect)
      {
         test = true;
      }
      cont = cont&&finish(test,"dump1");
   }
   if (cont)
   {
      start();
      t.select("test.tmp");
      t.dump(got,term);
      FakeData* tmp = dynamic_cast<FakeData*>(t.find("test.tmp"));
      bool test = tmp->touched;
      cont = cont&&finish(test,"dump2");
   }
   if (cont)
   {
      start();
      bool test = false;
      try
      {
         t.dump(got,term);
      }
      catch (int x)
      {
         test = x==1;
      }
      cont = cont&&finish(test,"dump3");
   }
   return cont;
}



bool unit::datamap::query()
{
   LinuxTerm term;
   GetOpts got("one");
   DataMap t;
   t.open("test.tmp","FakeData");
   bool cont = true;
   {
      start();
      bool test = false;
      try
      {
         t.query(got,term);
      }
      catch (DataMap::NoSelect)
      {
         test = true;
      }
      cont = cont&&finish(test,"query1");
   }
   if (cont)
   {
      start();
      t.select("test.tmp");
      t.query(got,term);
      FakeData* tmp = dynamic_cast<FakeData*>(t.find("test.tmp"));
      bool test = tmp->touched;
      cont = cont&&finish(test,"query2");
   }
   if (cont)
   {
      start();
      bool test = false;
      try
      {
         t.query(got,term);
      }
      catch (int x)
      {
         test = x==2;
      }
      cont = cont&&finish(test,"query3");
   }
   return cont;
}



bool unit::datamap::iterate()
{
   start();
   DataMap t;
   t.open("test1.tmp","FakeData");
   t.open("test2.tmp","FakeData");
   t.open("test3.tmp","FakeData");
   int count {0};
   for (auto i = t.begin();i!=t.end();++i)
   {
      ++count;
   }
   bool test = count==3;
   return finish(test,"iterate1");
}



bool unit::datamap::selected()
{
   bool cont {true};
   DataMap t;
   {
      start();
      bool test {t.selected()==t.end()};
      cont = cont&&finish(test,"selected1");
   }
   {
      start();
      t.open("test.tmp","FakeData",true);
      bool test {t.selected()==t.begin()};
      cont = cont&&finish(test,"selected2");
   }
   return cont;
}



bool unit::datamap::iter_file()
{
   start();
   std::string file("test.tmp");
   DataMap t;
   t.open(file,"FakeData");
   auto i = t.begin();
   bool test = i.file()==file;
   return finish(test,"iter_file1");
}



bool unit::datamap::iter_type()
{
   start();
   std::string type("FakeData");
   DataMap t;
   t.open("test.tmp",type);
   auto i = t.begin();
   bool test = i.type()==type;
   return finish(test,"iter_type1");
}
