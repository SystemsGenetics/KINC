#include "datamap.h"
#include "plugins/plugins.h"



DataPlugin& DataMap::open(const string& name, bool select)
{
   auto n = name.begin();
   while (n!=name.end()&&*n!=':')
   {
      ++n;
   }
   string file(name.begin(),n);
   string type(++n,name.end());
   bool cond = _map.find(file)==_map.end();
   assert<AlreadyExists>(cond,__FILE__,__LINE__);
   std::unique_ptr<DataPlugin> nd(KINCPlugins::new_data(type,file));
   assert<InvalidType>(bool(nd),__FILE__,__LINE__);
   auto x = _map.emplace(file,std::move(nd));
   auto i = x.first;
   if (select)
   {
      _i = i;
   }
   return *(i->second);
}



void DataMap::close(const string& file)
{
   auto i = get(file);
   if (_i==i)
   {
      _i = _map.end();
   }
   _map.erase(i);
}



void DataMap::load(GetOpts& ops, Terminal& tm)
{
   bool cond = _i!=_map.end();
   assert<NoSelect>(cond,__FILE__,__LINE__);
   try
   {
      _i->second->load(ops,tm);
   }
   catch (...)
   {
      _map.erase(_i);
      _i = _map.end();
      throw;
   }
}



void DataMap::dump(GetOpts& ops, Terminal& tm)
{
   bool cond = _i!=_map.end();
   assert<NoSelect>(cond,__FILE__,__LINE__);
   try
   {
      _i->second->dump(ops,tm);
   }
   catch (...)
   {
      _map.erase(_i);
      _i = _map.end();
      throw;
   }
}



void DataMap::query(GetOpts& ops, Terminal& tm)
{
   bool cond = _i!=_map.end();
   assert<NoSelect>(cond,__FILE__,__LINE__);
   try
   {
      _i->second->query(ops,tm);
   }
   catch (...)
   {
      _map.erase(_i);
      _i = _map.end();
      throw;
   }
}



DataMap::Map::iterator DataMap::get(const string& file)
{
   auto i = _map.find(file);
   bool cond = i!=_map.end();
   assert<DoesNotExist>(cond,__FILE__,__LINE__);
   return i;
}
