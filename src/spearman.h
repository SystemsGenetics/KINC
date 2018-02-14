#ifndef SPEARMAN_H
#define SPEARMAN_H
#include <ace/core/AceCore.h>
#include <ace/core/AceOpenCL.h>

#include "expressionmatrix.h"
#include "ccmatrix.h"
#include "correlationmatrix.h"
#include "genepair_vector.h"



class Spearman : public EAbstractAnalytic
{
   Q_OBJECT

public:
   ~Spearman();

   enum Arguments
   {
      InputData = 0
      ,ClusterData
      ,OutputData
      ,MinSamples
      ,MinExpression
      ,MinCorrelation
      ,MaxCorrelation
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
   virtual void runSerial() override final;
   virtual int getBlockSize() override;
   virtual bool runBlock(int block) override final;
   virtual void finish() override final {}

private:
   int fetchData(const GenePair::Vector& vector, const CCMatrix::Pair& pair, int k, double *a, double *b);

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

      Block(EOpenCLDevice& device, int size, int workSize, int kernelSize)
      {
         pairs = device.makeBuffer<cl_int2>(kernelSize).release();
         sampleMasks = device.makeBuffer<cl_char>(size * kernelSize).release();
         workBuffer = device.makeBuffer<cl_float>(2*workSize * kernelSize).release();
         rankBuffer = device.makeBuffer<cl_int>(2*workSize * kernelSize).release();
         results = device.makeBuffer<cl_float>(kernelSize).release();

         if ( !device )
         {
            E_MAKE_EXCEPTION(e);
            throw e;
         }
      }

      ~Block()
      {
         delete pairs;
         delete sampleMasks;
         delete workBuffer;
         delete rankBuffer;
         delete results;
      }

      int state {Start};
      GenePair::Vector vector;
      int cluster;
      EOpenCLEvent event;
      EOpenCLBuffer<cl_int2>* pairs;
      EOpenCLBuffer<cl_char>* sampleMasks;
      EOpenCLBuffer<cl_float>* workBuffer;
      EOpenCLBuffer<cl_int>* rankBuffer;
      EOpenCLBuffer<cl_float>* results;
   };

   void initializeKernel();
   void initializeBlockExpressions();
   void initializeKernelArguments();
   void runStartBlock(Block& block);
   void runLoadBlock(Block& block);
   void runExecuteBlock(Block& block);
   void runReadBlock(Block& block);

   ExpressionMatrix* _input {nullptr};
   CCMatrix* _cMatrix {nullptr};
   CorrelationMatrix* _output {nullptr};
   int _minSamples {30};
   float _minExpression {-INFINITY};
   float _minCorrelation {0.5};
   float _maxCorrelation {1.0};
   int _blockSize {4};
   int _kernelSize {4096};

   Block** _blocks {nullptr};
   EOpenCLProgram* _program {nullptr};
   EOpenCLKernel* _kernel {nullptr};
   EOpenCLBuffer<cl_float>* _expressions {nullptr};
   GenePair::Vector _vector;
   int _cluster {0};
   GenePair::Vector _nextVector;
   CCMatrix::Pair _inPair;
   CorrelationMatrix::Pair _outPair;
   qint64 _totalPairs;
   qint64 _pairsComplete {0};
   int _lastPercent {0};
};



#endif
