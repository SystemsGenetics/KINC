#ifndef EMATRIX_H
#define EMATRIX_H
#include "../../../dataplugin.h"



class ematrix : public DataPlugin
{
public:
   using string = std::string;
   ematrix(const string& type, const string& file): DataPlugin(type,file) {}
   void load(GetOpts &ops, Terminal &tm) override final {}
   void dump(GetOpts &ops, Terminal &tm) override final {}
   void query(GetOpts &ops, Terminal &tm) override final {}
   bool empty() override final { return true; }
};



#endif
