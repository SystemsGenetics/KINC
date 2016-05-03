#ifndef CMATRIX_HH
#define CMATRIX_HH
#include "../../../dataplugin.h"



class cmatrix : public DataPlugin
{
public:
   cmatrix(const string& type, const string& file): DataPlugin(type,file) {}
   void load(GetOpts &ops, Terminal &tm) override final {}
   void dump(GetOpts &ops, Terminal &tm) override final {}
   void query(GetOpts &ops, Terminal &tm) override final {}
   bool empty() override final { return true; }
};



#endif
