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
   virtual void setArgument(int argument, QVariant value) override final;
   virtual void setArgument(int argument, QFile* file) override final;
   virtual void setArgument(int argument, EAbstractData* data) override final;
   virtual quint32 getCapabilities() const override final { return Capabilities::Serial; }
   virtual bool initialize() override final;
   virtual void runSerial() override final;
   virtual int getBlockSize() override final { return 0; }
   virtual bool runBlock(int /*block*/) override final { return false; }
   virtual void finish() override final {}
private:
   float determineThreshold();
   float determineChi(float threshold);
   void generateGeneThresholds();
   float* generatePruneMatrix(float threshold, int* size);
   float* generateMatrixEigens(float* pruneMatrix, int* size);
   float getNNSDChiSquare(float* eigens, int size);
   CorrelationMatrix* _input;
   QFile* _output;
   float _initialThreshold;
   float _thresholdStep;
   float _thresholdMinimum;
};



#endif
