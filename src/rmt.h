#ifndef RMT_H
#define RMT_H
#include <ace.h>
#include "cmatrix.h"



namespace Ace = AccelCompEng;



class RMT : public Ace::Analytic
{
public:
   struct NoDataInput : public Ace::Exception { using Ace::Exception::Exception; };
   void input(Ace::Data*) override final {}
   void output(Ace::Data*) override final {}
protected:
   void execute_cl(Ace::GetOpts&,Ace::Terminal&) override final {}
   void execute_pn(Ace::GetOpts&,Ace::Terminal&) override final {}
private:
   double* unfolding(float* e, int size, int m);
   float* degenerate(float* eigens, int size, int* newSize);
   double getNNSDChiSquare(float* eigens, int size);
   double getNNSDPaceChiSquare(float* eigens, int size, double bin, int pace);
   static void swapD(double* l, int idx1, int idx2);
   static void quickSortD(double* l, int size);
   CMatrix* _in {nullptr};
   constexpr static int minEigenVectorSize {100};
   constexpr static float nnsdHistogramBin {0.05};
   constexpr static float chiSquareTestThreshold {99.607};
   constexpr static int minUnfoldingPace {10};
   constexpr static int maxUnfoldingPace {41};
};



#endif
