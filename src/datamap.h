#ifndef DATAMAP_H
#define DATAMAP_H
#include <map>
#include <string>
#include <memory>
#include "dataplugin.h"
#include "getopts.h"
#include "exception.h"



/// @brief Manages all data objects.
///
/// Manages list of all open data objects the program has loaded. Handles
/// passing user commands to individual data obejcts, along with selecting an
/// active data object, opening and closing data objects, and listing all loaded
/// data objects.
///
/// @warning There can only be one instance of this class in the program.
///
/// @author Josh Burns
/// @date 21 March 2016
class DataMap
{
public:
   // *
   // * EXCEPTIONS
   // *
   struct Exception;
   struct InvalidUse;
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
   DataMap();
   ~DataMap();
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
   DataPlugin* open(const string&,const string&,bool = false);
   bool close(const string&);
   void select(const string&);
   bool unselect();
   void load(GetOpts&,Terminal&);
   void dump(GetOpts&,Terminal&);
   void query(GetOpts&,Terminal&);
   DataPlugin* find(const string&);
   DataPlugin* current();
   Iterator begin();
   Iterator end();
   Iterator selected();
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
   // * STATIC VARIABLES
   // *
   /// Used to determine if an instance of this class exists or not.
   static bool _lock;
   // *
   // * VARIABLES
   // *
   /// Contains all loaded data objects.
   Map _map;
   /// Iterator that points to selected data object, or end of list iterator if
   /// no object is selected.
   Map::iterator _i;
};



/// Forward only iterator for listing all data objects in the DataMap class.
///
/// @author Josh Burns
/// @date 21 March 2016
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
   void operator++();
   bool operator!=(const Iterator&);
   bool operator==(const Iterator&);
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
   /// Points to current data object this iterator points to or end of list.
   Map::iterator _i;
};



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

struct DataMap::InvalidUse : public DataMap::Exception
{
   InvalidUse(const char* file, int line):
      Exception(file,line,"DataMap::InvalidUse")
   {}
};

struct DataMap::AlreadyExists : public DataMap::Exception
{
   AlreadyExists(const char* file, int line):
      Exception(file,line,"DataMap::AlreadyExists")
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
