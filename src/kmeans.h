#ifndef KMEANS_H
#define KMEANS_H
#include <ace/core/AceCore.h>
#include <ace/core/AceOpenCL.h>

#include "expressionmatrix.h"
#include "ccmatrix.h"
#include "genepair_kmeans.h"
#include "genepair_vector.h"



class KMeans : public EAbstractAnalytic
{
   Q_OBJECT

public:
   ~KMeans();

   enum Arguments
   {
      InputData = 0
      ,OutputData
      ,MinSamples
      ,MinClusters
      ,MaxClusters
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
   float computeBIC(const GenePair::KMeans& model, int N, int D);
   CCMatrix::Pair computePair(const QVector<GenePair::Vector2>& X);

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

      Block(EOpenCLDevice& device, int N, int K, int kernelSize)
      {
         pairs = device.makeBuffer<cl_int2>(1 * kernelSize).release();
         workX = device.makeBuffer<cl_float2>(N * kernelSize).release();
         workMu = device.makeBuffer<cl_float2>(K * kernelSize).release();
         worky1 = device.makeBuffer<cl_int>(N * kernelSize).release();
         worky2 = device.makeBuffer<cl_int>(N * kernelSize).release();
         resultKs = device.makeBuffer<cl_int>(1 * kernelSize).release();
         resultLabels = device.makeBuffer<cl_int>(N * kernelSize).release();
         if ( !device )
         {
            E_MAKE_EXCEPTION(e);
            throw e;
         }
      }

      ~Block()
      {
         delete pairs;
         delete workX;
         delete workMu;
         delete worky1;
         delete worky2;
         delete resultKs;
         delete resultLabels;
      }

      int state {Start};
      GenePair::Vector vector;
      EOpenCLEvent event;
      EOpenCLBuffer<cl_int2>* pairs;
      EOpenCLBuffer<cl_float2>* workX;
      EOpenCLBuffer<cl_float2>* workMu;
      EOpenCLBuffer<cl_int>* worky1;
      EOpenCLBuffer<cl_int>* worky2;
      EOpenCLBuffer<cl_int>* resultKs;
      EOpenCLBuffer<cl_int>* resultLabels;
   };

   void initializeKernel();
   void initializeBlockExpressions();
   void initializeKernelArguments();
   void runStartBlock(Block& block);
   void runLoadBlock(Block& block);
   void runExecuteBlock(Block& block);
   void runReadBlock(Block& block);

   ExpressionMatrix* _input {nullptr};
   CCMatrix* _output {nullptr};
   int _minSamples {30};
   int _minClusters {1};
   int _maxClusters {5};
   int _blockSize {4};
   int _kernelSize {4096};

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
