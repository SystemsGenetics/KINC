#ifndef DATA_H
#define DATA_H
#include <string>
#include "terminal.h"
#include "getopts.h"



class Data
{
public:
   virtual ~Data() = default;
   virtual void load(GetOpts&,Terminal&) = 0;
   virtual void dump(GetOpts&,Terminal&) = 0;
   virtual void query(GetOpts&,Terminal&) = 0;
   virtual bool empty() = 0;
};



#endif
