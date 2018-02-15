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
      ,MinExpression
      ,MinClusters
      ,MaxClusters
      ,NumInits
      ,MaxIterations
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
   void fetchData(const GenePair::Vector& vector, QVector<GenePair::Vector2>& X, QVector<cl_char>& labels);
   float computeBIC(int K, float logL, int N, int D);
   void computeModel(const QVector<GenePair::Vector2>& X, int& bestK, QVector<cl_char>& bestLabels);
   void savePair(const GenePair::Vector& vector, int K, const cl_char *labels, int N);

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
         work_X = device.makeBuffer<cl_float2>(N * kernelSize).release();
         work_labels = device.makeBuffer<cl_char>(3*N * kernelSize).release();
         work_means = device.makeBuffer<cl_float2>(K * kernelSize).release();
         result_K = device.makeBuffer<cl_int>(1 * kernelSize).release();
         result_labels = device.makeBuffer<cl_char>(N * kernelSize).release();

         if ( !device )
         {
            E_MAKE_EXCEPTION(e);
            throw e;
         }
      }

      ~Block()
      {
         delete pairs;
         delete work_X;
         delete work_labels;
         delete work_means;
         delete result_K;
         delete result_labels;
      }

      bool isWaiting()
      {
         for ( auto& event : events )
         {
            if ( !event.isDone() )
            {
               return true;
            }
         }

         return false;
      }

      void checkAllEvents()
      {
         for ( auto& event : events )
         {
            if ( !event )
            {
               E_MAKE_EXCEPTION(e);
               event.fillException(e);
               throw e;
            }
         }
      }

      int state {Start};
      GenePair::Vector vector;
      QVector<EOpenCLEvent> events;
      EOpenCLBuffer<cl_int2>* pairs;
      EOpenCLBuffer<cl_float2>* work_X;
      EOpenCLBuffer<cl_char>* work_labels;
      EOpenCLBuffer<cl_float2>* work_means;
      EOpenCLBuffer<cl_int>* result_K;
      EOpenCLBuffer<cl_char>* result_labels;
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
   float _minExpression {-INFINITY};
   int _minClusters {1};
   int _maxClusters {5};
   int _numInits {10};
   int _maxIterations {300};
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
