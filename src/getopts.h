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
   struct InvalidSyntax;
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
   const string& orig() const;
   int com_size() const;
   bool com_empty() const;
   int com_get(std::initializer_list<string>);
   string& com_front();
   void com_pop();
   int size() const;
   bool empty() const;
   bool has_opt(const string&,bool = false);
   Iterator begin();
   Iterator end();
   Iterator erase(Iterator);
private:
   // *
   // * DECLERATIONS
   // *
   using clist = std::list<string>;
   using oplist = std::list<std::pair<string,string>>;
   // *
   // * VARIABLES
   // *
   string _orig;
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
   const string& key() const;
   template<class T> T value() const;
   bool is_key(const string&) const;
   bool val_empty() const;
   // *
   // * OPERATORS
   // *
   void operator++();
   bool operator!=(const Iterator&);
   bool operator==(const Iterator&);
private:
   // *
   // * BASIC METHODS
   // *
   Iterator(const oplist::iterator&);
   // *
   // * VARIABLES
   // *
   oplist::iterator _i;
};



inline const GetOpts::string& GetOpts::orig() const
{
   return _orig;
}



inline int GetOpts::com_size() const
{
   return _comms.size();
}



inline bool GetOpts::com_empty() const
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



inline int GetOpts::size() const
{
   return _opts.size();
}



inline bool GetOpts::empty() const
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



inline const GetOpts::string& GetOpts::Iterator::key() const
{
   return _i->first;
}



template<class T> T GetOpts::Iterator::value() const
{
   T ret;
   using istring = std::istringstream;
   istring buffer(_i->second);
   bool cond = buffer >> ret;
   assert<InvalidType>(cond,__FILE__,__LINE__);
   return ret;
}



inline bool GetOpts::Iterator::is_key(const string& cmp) const
{
   return _i->first==cmp;
}



inline bool GetOpts::Iterator::val_empty() const
{
   return _i->second.empty();
}



inline void GetOpts::Iterator::operator++()
{
   ++_i;
}



inline bool GetOpts::Iterator::operator!=(const Iterator& cmp)
{
   return _i!=cmp._i;
}



inline bool GetOpts::Iterator::operator==(const Iterator& cmp)
{
   return _i==cmp._i;
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
      Exception(file,line,"GetOpts::InvalidType")
   {}
};

struct GetOpts::InvalidSyntax : public GetOpts::Exception
{
   InvalidSyntax(const char* file, int line):
      Exception(file,line,"GetOpts::InvalidSyntax")
   {}
};



#endif
