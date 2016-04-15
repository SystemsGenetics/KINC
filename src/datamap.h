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



/// Select a data object ot have focus.
///
/// @param file Name of loaded data object to select.
inline void DataMap::select(const string& file)
{
   _i = get(file);
}



/// Get beginning of list iterator.
///
/// @return Iterator.
inline DataMap::Iterator DataMap::begin()
{
   return _map.begin();
}



/// Get end of list iterator.
///
/// @return Iterator.
inline DataMap::Iterator DataMap::end()
{
   return _map.end();
}



/// Return iterator of currently selected data object, if any.
///
/// @return Iterator of currently selected object, end of list iterator if no
/// object is selected.
inline DataMap::Iterator DataMap::selected()
{
   return _i;
}



/// Get file name of iterator's data object.
///
/// @return File name of data object.
inline DataMap::Iterator::string DataMap::Iterator::file()
{
   return _i->first;
}



/// Get iterator's data object type.
///
/// @return Data object type.
inline DataMap::Iterator::string DataMap::Iterator::type()
{
   return _i->second->type();
}



/// Iterate to next data object in list of objects.
inline void DataMap::Iterator::operator++()
{
   ++_i;
}



/// Compare between two data object iterators.
///
/// @return False if they do not point to the same data object, else true.
inline bool DataMap::Iterator::operator!=(const Iterator& cmp)
{
   return _i!=cmp._i;
}



/// Compare between two data object iterators.
///
/// @return True if they point to the same data object, else false.
inline bool DataMap::Iterator::operator==(const Iterator& cmp)
{
   return _i==cmp._i;
}



/// Initializes data object iterator from internal container iterator.
///
/// @param i Internal container iterator.
inline DataMap::Iterator::Iterator(Map::iterator i):
   _i(i)
{}



/// Base exception class for any exception thrown from DataMap class.
struct DataMap::Exception : public ::Exception
{
   using ::Exception::Exception;
};

/// Exception thrown if a second instance of the DataMap class is created.
struct DataMap::InvalidUse : public DataMap::Exception
{
   InvalidUse(const char* file, int line):
      Exception(file,line,"DataMap::InvalidUse")
   {}
};

/// Exception thrown if a data object is opened with the same file name as
/// another data object already opened.
struct DataMap::AlreadyExists : public DataMap::Exception
{
   AlreadyExists(const char* file, int line):
      Exception(file,line,"DataMap::AlreadyExists")
   {}
};

/// Thrown if attempting to select a data object with a file name is not in
/// the list of opened data objects.
struct DataMap::DoesNotExist : public DataMap::Exception
{
   DoesNotExist(const char* file, int line):
      Exception(file,line,"DataMap::DoesNotExist")
   {}
};

/// Thrown if the commands load, dump, or query are called and no data object
/// is selected.
struct DataMap::NoSelect : public DataMap::Exception
{
   NoSelect(const char* file, int line):
      Exception(file,line,"DataMap::NoSelect")
   {}
};

/// Thrown is the data type given in the open command cannot be found using the
/// KINCPlugins system.
struct DataMap::InvalidType : public DataMap::Exception
{
   InvalidType(const char* file, int line):
      Exception(file,line,"DataMap::InvalidType")
   {}
};



#endif
