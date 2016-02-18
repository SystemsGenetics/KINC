#ifndef DATA_H
#define DATA_H
#include <string>



class Data
{
public:
   virtual ~Data() = default;
   virtual void option(const std::string&,const std::string&) = 0;
   virtual void load(const std::string&) = 0;
   virtual void dump(const std::string&) = 0;
   virtual void query() = 0;
};



#endif
