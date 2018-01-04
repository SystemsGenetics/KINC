#ifndef PEARSON_H
#define PEARSON_H
#include <ace/core/AceCore.h>
#include <ace/core/AceOpenCL.h>

#include "genepair_vector.h"



class ExpressionMatrix;
class CorrelationMatrix;



class Pearson : public EAbstractAnalytic
{
   Q_OBJECT
public:
   ~Pearson();
   enum Arguments
   {
      InputData = 0
      ,OutputData
      ,Minimum
      ,MinThreshold
      ,MaxThreshold
      ,BlockSize
      ,KernelSize
      ,Total
   };
   virtual int getArgumentCount() override final { return Total; }
   virtual ArgumentType getArgumentData(int argument) override final;
   virtual QVariant getArgumentData(int argument, Role role) override final;
   virtual void setArgument(int argument, QVariant value) override final;
   virtual void setArgument(int argument, EAbstractData* data) override final;
   quint32 getCapabilities() const override final
      { return Capabilities::Serial|Capabilities::OpenCL; }
   virtual bool initialize() override final;
   virtual int getBlockSize() override;
   virtual bool runBlock(int block) override final;
   virtual void finish() override final {}
   virtual void runSerial() override final;
private:
   struct Block
   {
      enum
      {
         Start
         ,Load
         ,Execute
         ,Read
         ,Done
      };
      Block(EOpenCLDevice& device, int size, int kernelSize)
      {
         references = device.makeBuffer<cl_int>(2*kernelSize).release();
         answers = device.makeBuffer<cl_float>(kernelSize).release();
         workBuffer = device.makeBuffer<cl_float>(2*size*kernelSize).release();
         if ( !device )
         {
            E_MAKE_EXCEPTION(e);
            throw e;
         }
      }
      ~Block()
      {
         delete references;
         delete answers;
         delete workBuffer;
      }
      int state {Start};
      GenePair::Vector vector;
      EOpenCLEvent event;
      EOpenCLBuffer<cl_int>* references;
      EOpenCLBuffer<cl_float>* answers;
      EOpenCLBuffer<cl_float>* workBuffer;
   };
   void initializeKernel();
   void initializeBlockExpressions();
   void initializeKernelArguments();
   void runStartBlock(Block& block);
   void runLoadBlock(Block& block);
   void runExecuteBlock(Block& block);
   void runReadBlock(Block& block);
   ExpressionMatrix* _input {nullptr};
   CorrelationMatrix* _output {nullptr};
   int _minimum {30};
   int _blockSize {4};
   int _kernelSize {4096};
   float _minThreshold {0.5};
   float _maxThreshold {1.0};
   Block** _blocks {nullptr};
   EOpenCLProgram* _program {nullptr};
   EOpenCLKernel* _kernel {nullptr};
   EOpenCLBuffer<cl_float>* _expressions {nullptr};
   GenePair::Vector _vector;
   GenePair::Vector _nextVector;
   qint64 _totalPairs;
   qint64 _pairsComplete {0};
   int _lastPercent {0};
};



#endif