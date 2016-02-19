#ifndef DATA_H
#define DATA_H
#include <string>
#include "terminal.h"



class Data
{
public:
   virtual ~Data() = default;
   virtual void option(const std::string&,const std::string&) = 0;
   virtual void load(const std::string&,Terminal&) = 0;
   virtual void dump(const std::string&,Terminal&) = 0;
   virtual void query(Terminal&) = 0;
};



#endif
