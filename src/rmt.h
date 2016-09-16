#ifndef RMT_H
#define RMT_H
#include <ace.h>
#include "cmatrix.h"



namespace Ace = AccelCompEng;



class RMT : public Ace::Analytic
{
public:
   struct NoDataInput : public Ace::Exception { using Ace::Exception::Exception; };
   void input(Ace::Data*) override final;
   void output(Ace::Data*) override final;
protected:
   void execute_cl(Ace::GetOpts&,Ace::Terminal&) override final;
   void execute_pn(Ace::GetOpts&,Ace::Terminal&) override final;
private:
   CMatrix* _in {nullptr};
};



#endif
