#include "datamap.h"
#include "unit.h"
#include "plugins/plugins.h"



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
   try
   {
      ret = construct()&&
            open()&&
            close();
   }
   catch (...)
   {
      end();
      throw;
   }
   end();
   return ret;
}



bool unit::datamap::construct()
{
   start();
   DataMap t;
   return finish(true,"construct1");
}



bool unit::datamap::open()
{
   DataMap t;
   bool cont = true;
   {
      start();
      DataPlugin* n = t.open("test:FakeData",true);
      bool test = n==t.find("test");
      cont = cont&&finish(test,"open1");
   }
   if (cont)
   {
      start();
      bool test = false;
      try
      {
         t.open("test2:UnknownData");
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
         t.open("test:FakeData");
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
   start();
   DataMap t;
   t.open("test:FakeData");
   t.close("test");
   bool test = false;
   try
   {
      t.find("test");
   }
   catch (DataMap::DoesNotExist)
   {
      test = true;
   }
   return finish(test,"close1");
}
