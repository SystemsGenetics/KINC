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
   struct Exception;
   struct InvalidType;
   class Iterator;
   using string = std::string;
   using clist = std::list<string>;
   using oplist = std::list<std::pair<string,string>>;
   GetOpts(const string&);
   int get_com(std::initializer_list<string>);
   void pop_com();
   int size();
   bool empty();
   Iterator begin();
   Iterator end();
   Iterator erase(Iterator);
private:
   inline clist explode(const string&);
   clist _comms;
   oplist _opts;
};



class GetOpts::Iterator
{
public:
   friend class GetOpts;
   using string = GetOpts::string;
   using oplist = GetOpts::oplist;
   const string& key();
   Iterator& operator>>(short&);
   Iterator& operator>>(unsigned short&);
   Iterator& operator>>(int&);
   Iterator& operator>>(unsigned int&);
   Iterator& operator>>(long&);
   Iterator& operator>>(unsigned long&);
   Iterator& operator>>(float&);
   Iterator& operator>>(double&);
   Iterator& operator>>(string&);
   void operator++();
   bool operator!=(const Iterator&);
private:
   using istring = std::istringstream;
   Iterator(const oplist::iterator&);
   oplist::iterator _i;
};



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
