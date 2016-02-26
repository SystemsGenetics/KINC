#ifndef DATAMAP_H
#define DATAMAP_H
#include <unordered_map>
#include <string>
#include <memory>
#include "dataplugin.h"
#include "getopts.h"
#include "exception.h"



class DataMap
{
public:
   // *
   // * EXCEPTIONS
   // *
   struct Exception;
   struct AlreadyExists;
   struct DoesNotExist;
   struct NoSelect;
   struct InvalidType;
   // *
   // * DECLERATIONS
   // *
   class Iterator;
   using string = std::string;
   using wptr = std::weak_ptr<DataPlugin>;
   // *
   // * BASIC METHODS
   // *
   DataMap(const DataMap&) = delete;
   DataMap(DataMap&&) = delete;
   DataMap& operator=(const DataMap&) = delete;
   DataMap& operator=(DataMap&&) = delete;
   inline DataMap();
   // *
   // * FUNCTIONS
   // *
   wptr open(const string&,bool = false);
   void close(const string&);
   inline void select(const string&);
   void load(GetOpts&,Terminal&);
   void dump(GetOpts&,Terminal&);
   void query(GetOpts&,Terminal&);
   inline wptr find(const string&);
   inline Iterator begin();
   inline Iterator end();
private:
   // *
   // * DECLERATIONS
   // *
   using Map = std::unordered_map<std::string,std::shared_ptr<DataPlugin>>;
   // *
   // * FUNCTIONS
   // *
   Map::iterator get(const string&);
   // *
   // * VARIABLES
   // *
   Map _map;
   Map::iterator _i;
};



class DataMap::Iterator
{
public:
   // *
   // * DECLERATIONS
   // *
   friend class DataMap;
   using string = DataMap::string;
   // *
   // * FUNCTIONS
   // *
   inline string file();
   inline string type();
   // *
   // * OPERATORS
   // *
   inline void operator++();
   inline bool operator!=(const Iterator&);
private:
   // *
   // * DECLERATIONS
   // *
   using Map = DataMap::Map;
   // *
   // * BASIC METHODS
   // *
   inline Iterator(Map::iterator);
   // *
   // * VARIABLES
   // *
   Map::iterator _i;
};



inline DataMap::DataMap():
   _i {_map.end()}
{}



inline void DataMap::select(const string& file)
{
   _i = get(file);
}


inline DataMap::wptr DataMap::find(const string& file)
{
   return get(file)->second;
}



inline DataMap::Iterator DataMap::begin()
{
   return _map.begin();
}



inline DataMap::Iterator DataMap::end()
{
   return _map.end();
}



inline DataMap::Iterator::string DataMap::Iterator::file()
{
   return _i->first;
}



inline DataMap::Iterator::string DataMap::Iterator::type()
{
   return _i->second->type();
}



inline void DataMap::Iterator::operator++()
{
   ++_i;
}



inline bool DataMap::Iterator::operator!=(const Iterator& cmp)
{
   return _i!=cmp._i;
}



inline DataMap::Iterator::Iterator(Map::iterator i):
   _i(i)
{}



struct DataMap::Exception : public ::Exception
{
   using ::Exception::Exception;
};

struct DataMap::AlreadyExists : public DataMap::Exception
{
   AlreadyExists(const char* file, int line):
      Exception(file,line,"DataMap::InvalidUse")
   {}
};

struct DataMap::DoesNotExist : public DataMap::Exception
{
   DoesNotExist(const char* file, int line):
      Exception(file,line,"DataMap::DoesNotExist")
   {}
};

struct DataMap::NoSelect : public DataMap::Exception
{
   NoSelect(const char* file, int line):
      Exception(file,line,"DataMap::NoSelect")
   {}
};

struct DataMap::InvalidType : public DataMap::Exception
{
   InvalidType(const char* file, int line):
      Exception(file,line,"DataMap::InvalidType")
   {}
};



#endif
