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
   int com_size();
   bool com_empty();
   int com_get(std::initializer_list<string>);
   string& com_front();//NOT TESTED!!
   void com_pop();
   int size();
   bool empty();
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
   const string& key();
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
   void operator++();
   bool operator==(const Iterator&);
   bool operator!=(const Iterator&);
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
