#ifndef RMT_H
#define RMT_H
#include <ace.h>
#include "cmatrix.h"



namespace Ace = AccelCompEng;



class RMT : public Ace::Analytic
{
public:
   struct InvalidDataType : public Ace::Exception { using Ace::Exception::Exception; };
   struct TooManyInputs : public Ace::Exception { using Ace::Exception::Exception; };
   struct TooManyOutputs : public Ace::Exception { using Ace::Exception::Exception; };
   struct NoDataInput : public Ace::Exception { using Ace::Exception::Exception; };
   void input(Ace::Data*) override final;
   void output(Ace::Data*) override final;
protected:
   void execute_cl(Ace::GetOpts&,Ace::Terminal&) override final;
   void execute_pn(Ace::GetOpts&,Ace::Terminal&) override final;
private:
   std::unique_ptr<float> matrix_eigens(const std::unique_ptr<double>& matrix, int size);
   std::unique_ptr<double> prune_matrix(float threshold, int& size);
   void generate_gene_thresholds(Ace::Terminal&);
   double* unfolding(float* e, int size, int m);
   float* degenerate(float* eigens, int size, int* newSize);
   double getNNSDChiSquare(float* eigens, int size);
   double getNNSDPaceChiSquare(float* eigens, int size, double bin, int pace);
   static void swapD(double* l, int idx1, int idx2);
   static void quickSortD(double* l, int size);
   CMatrix* _in {nullptr};
   constexpr static int minEigenVectorSize {100};
   constexpr static float nnsdHistogramBin {0.05};
   constexpr static float _chiSquareTestThreshold {0.99607};
   constexpr static float _chiMinimum {0.5};
   constexpr static float _chiStep {0.001};
   constexpr static int minUnfoldingPace {10};
   constexpr static int maxUnfoldingPace {41};
   std::unique_ptr<float> _geneThresh;
};



#endif
