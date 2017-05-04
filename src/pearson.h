#ifndef PEARSON_H
#define PEARSON_H
#include <ace.h>
#include "ematrix.h"
#include "cmatrix.h"



namespace Ace = AccelCompEng;



class Pearson : public Ace::Analytic
{
public:
   struct InvalidDataType : public Ace::Exception { using Ace::Exception::Exception; };
   struct TooManyInputs : public Ace::Exception { using Ace::Exception::Exception; };
   struct TooManyOutputs : public Ace::Exception { using Ace::Exception::Exception; };
   struct NoDataInput : public Ace::Exception { using Ace::Exception::Exception; };
   struct NoDataOutput : public Ace::Exception { using Ace::Exception::Exception; };
   struct TooManySamples : public Ace::Exception { using Ace::Exception::Exception; };
   void input(Ace::Data*) override final;
   void output(Ace::Data*) override final;
protected:
   void execute_cl(Ace::GetOpts&,Ace::Terminal&) override final;
   void execute_pn(Ace::GetOpts&,Ace::Terminal&) override final;
private:
   using elist = Ace::CLBuffer<cl_float>;
   int pow2_ceil(int);
   int pow2_floor(int);
   void calculate(Ace::Terminal& tm, Ace::CLKernel& kern, elist& expList, int size, int wSize,
                  int chunk, int blSize, int smSize, int minSize);
   EMatrix* _in {nullptr};
   CMatrix* _out {nullptr};
};



class PearsonHelpItem : public AbstractHelpItem
{
public:
   std::string getName() const override final
   {
      return "pearson";
   }
   std::string getDescription() const override final
   {
      return "pearson Analytic\n"
            "\n"
            "This analytic conducts pearson analysis on emx data and outputs cmx data.\n"
             "\n"
             "pearson --in=[inname] --out=[outname]:cmx (--slots|--bsize|--minsize)\n"
             "\n"
             "Takes data object [inname] that must be emx and outputs new file [outname] that is"
             " a cmx data object. --slots specifies how many parallel slots are processed at the"
             " same time. --bsize specifies how many kernels are ran in parallel per slot."
             " --minsize specifies he minimum number of samples a gene pair must share to be given"
             " a pearson coefficient.\n";
   }
};


#endif
