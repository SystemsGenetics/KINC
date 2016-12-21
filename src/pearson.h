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



#endif
