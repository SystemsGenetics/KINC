#ifndef DATAMAP_H
#define DATAMAP_H
#include <map>
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
   // *
   // * BASIC METHODS
   // *
   inline DataMap();//
   // *
   // * COPY METHODS
   // *
   DataMap(const DataMap&) = delete;
   DataMap& operator=(const DataMap&) = delete;
   // *
   // * MOVE METHODS
   // *
   DataMap(DataMap&&) = delete;
   DataMap& operator=(DataMap&&) = delete;
   // *
   // * FUNCTIONS
   // *
   DataPlugin* open(const string&,const string&,bool = false);//
   bool close(const string&);// RETEST NEEDED!!!
   void select(const string&);//
   bool unselect();// NOT TESTED!!!
   void load(GetOpts&,Terminal&);//
   void dump(GetOpts&,Terminal&);//
   void query(GetOpts&,Terminal&);//
   DataPlugin* find(const string&);//
   DataPlugin* current();//NOT TESTED!!!
   Iterator begin();//
   Iterator end();//
   Iterator selected();//NOT TESTED!!!
private:
   // *
   // * DECLERATIONS
   // *
   using Map = std::map<std::string,std::unique_ptr<DataPlugin>>;
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
   string file();
   string type();
   // *
   // * OPERATORS
   // *
   void operator++();//
   bool operator!=(const Iterator&);//
   bool operator==(const Iterator&);//NOT TESTED!!!
private:
   // *
   // * DECLERATIONS
   // *
   using Map = DataMap::Map;
   // *
   // * BASIC METHODS
   // *
   Iterator(Map::iterator);
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



inline DataMap::Iterator DataMap::begin()
{
   return _map.begin();
}



inline DataMap::Iterator DataMap::end()
{
   return _map.end();
}



inline DataMap::Iterator DataMap::selected()
{
   return _i;
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



inline bool DataMap::Iterator::operator==(const Iterator& cmp)
{
   return _i==cmp._i;
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
