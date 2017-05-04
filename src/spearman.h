#ifndef SPEARMAN_H
#define SPEARMAN_H
#include <ace.h>
#include "ematrix.h"
#include "cmatrix.h"



namespace Ace = AccelCompEng;



/// Spearman plugin class that uses spearman analysis to generate correlations between genes.
class Spearman : public Ace::Analytic
{
public:
#ifdef UNIT_TEST
   static void runUnitTests(Ace::Console& console);
#endif
   struct InvalidDataType : public Ace::Exception { using Ace::Exception::Exception; };
   struct TooManyInputs : public Ace::Exception { using Ace::Exception::Exception; };
   struct TooManyOutputs : public Ace::Exception { using Ace::Exception::Exception; };
   struct NoDataInput : public Ace::Exception { using Ace::Exception::Exception; };
   struct NoDataOutput : public Ace::Exception { using Ace::Exception::Exception; };
   struct TooManySamples : public Ace::Exception { using Ace::Exception::Exception; };
   /// ACE load input data object call.
   void input(Ace::Data*) override final;
   /// ACE load output data object call.
   void output(Ace::Data*) override final;
protected:
   /// ACE OpenCL accelerated execution call.
   void execute_cl(Ace::GetOpts&,Ace::Terminal&) override final;
   /// ACE CPU non-accelerated execution call.
   void execute_pn(Ace::GetOpts&,Ace::Terminal&) override final;
private:
   void theMaddening();
   void theMaddening2();
   using elist = Ace::CLBuffer<cl_float>;
   int pow2_ceil(int);
   int pow2_floor(int);
   void calculate(Ace::Terminal& tm, Ace::CLKernel& kern, elist& expList, int size, int blSize
                  , int smSize, int minSize);
   EMatrix* _in {nullptr};
   CMatrix* _out {nullptr};
};



class SpearmanHelpItem : public AbstractHelpItem
{
public:
   std::string getName() const override final
   {
      return "spearman";
   }
   std::string getDescription() const override final
   {
      return "spearman Analytic\n"
            "\n"
            "This analytic conducts spearman analysis on emx data and outputs cmx data.\n"
             "\n"
             "spearman --in=[inname] --out=[outname]:cmx (--slots|--bsize|--minsize)\n"
             "\n"
             "Takes data object [inname] that must be emx and outputs new file [outname] that is"
             " a cmx data object. --slots specifies how many parallel slots are processed at the"
             " same time. --bsize specifies how many kernels are ran in parallel per slot."
             " --minsize specifies he minimum number of samples a gene pair must share to be given"
             " a spearman coefficient.\n";
   }
};


#endif
