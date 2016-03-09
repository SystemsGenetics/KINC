#include "datamap.h"
#include "plugins/plugins.h"



DataPlugin* DataMap::open(const string& file, const string& type, bool select)
{
   using uptr = std::unique_ptr<DataPlugin>;
   bool cond = _map.find(file)==_map.end();
   assert<AlreadyExists>(cond,__FILE__,__LINE__);
   uptr nd(KINCPlugins::new_data(type,file));
   assert<InvalidType>(bool(nd),__FILE__,__LINE__);
   auto x = _map.emplace(file,std::move(nd));
   auto i = x.first;
   if (select)
   {
      _i = i;
   }
   return i->second.get();
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



DataPlugin* DataMap::find(const string& file)
{
   try
   {
      return get(file)->second.get();
   }
   catch (DoesNotExist)
   {
      return nullptr;
   }
}



DataMap::Map::iterator DataMap::get(const string& file)
{
   auto i = _map.find(file);
   bool cond = i!=_map.end();
   assert<DoesNotExist>(cond,__FILE__,__LINE__);
   return i;
}
