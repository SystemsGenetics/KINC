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
   void load(GetOpts&,Terminal&) override final {}
   void dump(GetOpts&,Terminal&) override final {}
   void query(GetOpts&,Terminal&) override final {}
   bool empty() override final { return true; }
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
      ret = true;
   }
   catch (...)
   {
      end();
      throw;
   }
   end();
   return ret;
}
