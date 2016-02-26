#ifndef DATAMAP_H
#define DATAMAP_H
#include <unordered_map>
#include <string>
#include <memory>
#include "dataplugin.h"



class DataMap
{
public:
   using smap = std::unordered_map<std::string,DataPlugin*>;
   // *
   // * BASIC METHODS
   // *
   DataMap(const DataMap&) = delete;
   DataMap(DataMap&&) = delete;
   DataMap& operator=(const DataMap&) = delete;
   DataMap& operator=(DataMap&&) = delete;
   DataMap() = default;
   // *
   // * FUNCTIONS
   // *
   bool add(const std::string&,DataPlugin*);
   bool del(const std::string&);
   bool del(DataPlugin*);
   bool exist(const std::string&);
   DataPlugin* find(const std::string&);
   smap::iterator begin();
   smap::iterator end();
private:
   // *
   // * VARIABLES
   // *
   smap _map;
};



#endif
