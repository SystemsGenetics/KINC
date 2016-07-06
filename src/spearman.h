#ifndef SPEARMAN_H
#define SPEARMAN_H
#include <ace.h>
#include "ematrix.h"
#include "cmatrix.h"



class spearman : public AccelCompEng::AnalyticPlugin
{
public:
   ACE_EXCEPTION(spearman,InvalidInputType)
   ACE_EXCEPTION(spearman,TooManyOutputs)
   ACE_EXCEPTION(spearman,NotEnoughInputs)
   ACE_EXCEPTION(spearman,NoOutput)
   ACE_EXCEPTION(spearman,DifferentSampleSizes)
   using GetOpts = AccelCompEng::GetOpts;
   using Terminal = AccelCompEng::Terminal;
   using FileMem = AccelCompEng::FileMem;
   using FileSync = AccelCompEng::FileSync;
   using KincFile = AccelCompEng::File;
   using FString = AccelCompEng::FString;
   using DataPlugin = AccelCompEng::DataPlugin;
   void input(DataPlugin*) override final;
   void output(DataPlugin*) override final;
protected:
   void execute_cl(GetOpts&,Terminal&) override final;
   void execute_pn(GetOpts&,Terminal&) override final;
private:
   std::vector<EMatrix*> _in;
   CMatrix* _out {nullptr};
};



#endif
