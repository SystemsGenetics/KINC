#ifndef GETOPTS_H
#define GETOPTS_H
#include <string>
#include <sstream>
#include <list>
#include <initializer_list>
#include "exception.h"



class GetOpts
{
public:
   // *
   // * EXCEPTIONS
   // *
   struct Exception;
   struct InvalidType;
   // *
   // * DECLERATIONS
   // *
   class Iterator;
   using string = std::string;
   // *
   // * BASIC METHODS
   // *
   GetOpts(const string&);
   // *
   // * FUNCTIONS
   // *
   inline int com_size();
   inline bool com_empty();
   int com_get(std::initializer_list<string>);
   inline string& com_front();//NOT TESTED!!
   inline void com_pop();
   inline int size();
   inline bool empty();
   inline Iterator begin();
   inline Iterator end();
   inline Iterator erase(Iterator);
private:
   // *
   // * DECLERATIONS
   // *
   using clist = std::list<string>;
   using oplist = std::list<std::pair<string,string>>;
   // *
   // * FUNCTIONS
   // *
   inline clist explode(const string&);
   // *
   // * VARIABLES
   // *
   clist _comms;
   oplist _opts;
};



class GetOpts::Iterator
{
public:
   // *
   // * DECLERATIONS
   // *
   friend class GetOpts;
   using string = GetOpts::string;
   using oplist = GetOpts::oplist;
   // *
   // * FUNCTIONS
   // *
   inline const string& key();
   inline bool empty();// NOT TESTED
   // *
   // * OPERATORS
   // *
   Iterator& operator>>(short&);
   Iterator& operator>>(unsigned short&);
   Iterator& operator>>(int&);
   Iterator& operator>>(unsigned int&);
   Iterator& operator>>(long&);
   Iterator& operator>>(unsigned long&);
   Iterator& operator>>(float&);
   Iterator& operator>>(double&);
   Iterator& operator>>(string&);
   inline void operator++();
   inline bool operator==(const Iterator&);
   inline bool operator!=(const Iterator&);
private:
   // *
   // * DECLERATIONS
   // *
   using istring = std::istringstream;
   // *
   // * BASIC METHODS
   // *
   inline Iterator(const oplist::iterator&);
   // *
   // * VARIABLES
   // *
   oplist::iterator _i;
};



inline int GetOpts::com_size()
{
   return _comms.size();
}



inline bool GetOpts::com_empty()
{
   return _comms.empty();
}



inline GetOpts::string& GetOpts::com_front()
{
   return _comms.front();
}



inline void GetOpts::com_pop()
{
   if (!_comms.empty())
   {
      _comms.pop_front();
   }
}



inline int GetOpts::size()
{
   return _opts.size();
}



inline bool GetOpts::empty()
{
   return _opts.empty();
}



inline GetOpts::Iterator GetOpts::begin()
{
   return {_opts.begin()};
}



inline GetOpts::Iterator GetOpts::end()
{
   return {_opts.end()};
}



inline GetOpts::Iterator GetOpts::erase(Iterator i)
{
   return {_opts.erase(i._i)};
}



inline const GetOpts::string& GetOpts::Iterator::key()
{
   return _i->first;
}



inline bool GetOpts::Iterator::empty()
{
   return _i->second.empty();
}



inline void GetOpts::Iterator::operator++()
{
   ++_i;
}



inline bool GetOpts::Iterator::operator==(const Iterator& cmp)
{
   return _i==cmp._i;
}



inline bool GetOpts::Iterator::operator!=(const Iterator& cmp)
{
   return _i!=cmp._i;
}



inline GetOpts::Iterator::Iterator(const oplist::iterator& i):
   _i(i)
{}



struct GetOpts::Exception : public ::Exception
{
   using ::Exception::Exception;
};

struct GetOpts::InvalidType : public GetOpts::Exception
{
   InvalidType(const char* file, int line):
      Exception(file,line,"Terminal::InvalidType")
   {}
};



#endif
