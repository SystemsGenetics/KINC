#ifndef RMT_H
#define RMT_H
#include <ace/core/AceCore.h>



class CorrelationMatrix;



class RMT : public EAbstractAnalytic
{
public:
   enum Arguments
   {
      InputData = 0
      ,OutputFile
      ,Total
   };
   virtual int getArgumentCount() override final { return Total; }
   virtual ArgumentType getArgumentData(int argument) override final;
   virtual QVariant getArgumentData(int argument, Role role) override final;
   virtual void setArgument(int argument, QFile* file) override final;
   virtual void setArgument(int argument, EAbstractData* data) override final;
   virtual quint32 getCapabilities() const override final { return Capabilities::Serial; }
   virtual bool initialize() override final;
   virtual void runSerial() override final;
private:
   float determineThreshold();
   float determineChi(float threshold);
   void generateGeneThresholds();
   double* generatePruneMatrix(float threshold, int* size);
   static float* generateMatrixEigens(double* pruneMatrix, int size);
   // FOR ALL FUNCTIONS BELOW
   // these were all pulled from KINC version 1
   static void swapF(float* l, int idx1, int idx2);
   static void quickSortF(float* l, int size);
   static void swapD(double* l, int idx1, int idx2);
   static void quickSortD(double* l, int size);
   static double* unfolding(float* e, int size, int m);
   static float* degenerate(float* eigens, int size, int* newSize);
   static double getNNSDChiSquare(float* eigens, int size);
   static double getNNSDPaceChiSquare(float* eigens, int size, double bin, int pace);
   // END OF ALL FUNCTIONS
   // FOR ALL STATIC VARIABLES BELOW
   // these were all pulled fron KINC version 1 used by its functions
   constexpr static float nnsdHistogramBin {0.05};
   constexpr static int minUnfoldingPace {10};
   constexpr static int maxUnfoldingPace {41};
   constexpr static int minEigenVectorSize {100};
   // END OF ALL STATIC VARIABLES
   CorrelationMatrix* _input {nullptr};
   QFile* _output {nullptr};
   float _initialThreshold {0.99607};
   float _thresholdStep {0.001};
   float _thresholdMinimum {0.5};
   std::unique_ptr<float> _geneThresholds;
};



#endif
