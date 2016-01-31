#ifndef DATA_H
#define DATA_H
#include <string>



class Data
{
public:
   virtual ~Data() = default;
   virtual void option(const std::string&,const std::string&) = 0;
   virtual bool load(const std::string&) = 0;
   virtual bool dump(const std::string&) = 0;
   virtual bool query() = 0;
   virtual void flush() = 0;
};



#endif
