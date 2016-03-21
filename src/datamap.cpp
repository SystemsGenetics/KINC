#include "datamap.h"
#include "plugins/plugins.h"



bool DataMap::_lock {false};



DataMap::DataMap():
   _i {_map.end()}
{
   assert<InvalidUse>(!_lock,__FILE__,__LINE__);
   _lock = true;
}



DataMap::~DataMap()
{
   _lock = false;
}



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



bool DataMap::close(const string& file)
{
   bool ret = false;
   auto i = get(file);
   if (_i==i)
   {
      ret = true;
      _i = _map.end();
   }
   _map.erase(i);
   return ret;
}



bool DataMap::unselect()
{
   bool ret = false;
   if (_i!=_map.end())
   {
      ret = true;
      _i = _map.end();
   }
   return ret;
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



DataPlugin* DataMap::current()
{
   bool cond = _i!=_map.end();
   assert<NoSelect>(cond,__FILE__,__LINE__);
   return _i->second.get();
}



DataMap::Map::iterator DataMap::get(const string& file)
{
   auto i = _map.find(file);
   bool cond = i!=_map.end();
   assert<DoesNotExist>(cond,__FILE__,__LINE__);
   return i;
}
