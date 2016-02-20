#include "datamap.h"



bool DataMap::add(const std::string& name, DataPlugin* data)
{
   bool ret = false;
   if (_map.find(name)==_map.end())
   {
      _map[name] = data;
      ret = true;
   }
   return ret;
}



bool DataMap::del(const std::string& name)
{
   bool ret = false;
   auto i = _map.find(name);
   if (i!=_map.end())
   {
      delete i->second;
      _map.erase(i);
      ret = true;
   }
   return ret;
}



bool DataMap::del(DataPlugin* data)
{
   bool ret = false;
   for (auto i = _map.begin();i!=_map.end();)
   {
      if (i->second==data)
      {
         delete i->second;
         _map.erase(i);
         ret = true;
      }
      else
      {
         ++i;
      }
   }
   return ret;
}



bool DataMap::exist(const std::string& name)
{
   return _map.find(name)!=_map.end();
}



DataPlugin* DataMap::find(const std::string& name)
{
   DataPlugin* ret {nullptr};
   auto i = _map.find(name);
   if (i!=_map.end())
   {
      ret = i->second;
   }
   return ret;
}



DataMap::smap::iterator DataMap::begin()
{
   return _map.begin();
}



DataMap::smap::iterator DataMap::end()
{
   return _map.end();
}
