#ifndef SPEARMAN_H
#define SPEARMAN_H
#include <ace.h>
#include "ematrix.h"
#include "cmatrix.h"



class Spearman : public AccelCompEng::AnalyticPlugin
{
public:
   ACE_EXCEPTION(Spearman,InvalidDataType)
   ACE_EXCEPTION(Spearman,TooManyInputs)
   ACE_EXCEPTION(Spearman,TooManyOutputs)
   ACE_EXCEPTION(Spearman,NoDataInput)
   ACE_EXCEPTION(Spearman,NoDataOutput)
   ACE_EXCEPTION(Spearman,TooManySamples)
   //ACE_EXCEPTION(spearman,NotEnoughInputs)
   //ACE_EXCEPTION(spearman,NoOutput)
   //ACE_EXCEPTION(spearman,DifferentSampleSizes)
   using GetOpts = AccelCompEng::GetOpts;
   using Terminal = AccelCompEng::Terminal;
   using FileMem = AccelCompEng::FileMem;
   using FSync = AccelCompEng::FileSync;
   using File = AccelCompEng::File;
   using FString = AccelCompEng::FString;
   using DataPlugin = AccelCompEng::DataPlugin;
   using CLContext = AccelCompEng::CLContext;
   using CLCommandQueue = AccelCompEng::CLCommandQueue;
   using CLEvent = AccelCompEng::CLEvent;
   using CLProgram = AccelCompEng::CLProgram;
   using CLKernel = AccelCompEng::CLKernel;
   void input(DataPlugin*) override final;
   void output(DataPlugin*) override final;
protected:
   void execute_cl(GetOpts&,Terminal&) override final;
   void execute_pn(GetOpts&,Terminal&) override final;
private:
   using elist = AccelCompEng::CLBuffer<cl_float>;
   int pow2_ceil(int);
   int pow2_floor(int);
   void calculate(Terminal&,CLKernel&,elist&,int,int,int);
   EMatrix* _in {nullptr};
   CMatrix* _out {nullptr};
};



#endif
