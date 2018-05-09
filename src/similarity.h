#ifndef SIMILARITY_H
#define SIMILARITY_H
#include <ace/core/AceCore.h>
#include <ace/core/AceOpenCL.h>

#include "expressionmatrix.h"
#include "ccmatrix.h"
#include "correlationmatrix.h"
#include "pairwise_clustering.h"
#include "pairwise_correlation.h"
#include "pairwise_gmm.h"



class Similarity : public EAbstractAnalytic
{
   Q_OBJECT

public:
   enum Arguments
   {
      InputData = 0
      ,ClusterData
      ,CorrelationData
      ,ClusteringArg
      ,CorrelationArg
      ,MinExpression
      ,MinSamples
      ,MinClusters
      ,MaxClusters
      ,CriterionArg
      ,RemovePreOutliers
      ,RemovePostOutliers
      ,MinCorrelation
      ,MaxCorrelation
      ,BlockSize
      ,KernelSize
      ,Total
   };

   enum class ClusteringMethod
   {
      None
      ,GMM
      ,KMeans
   };

   enum class CorrelationMethod
   {
      Pearson
      ,Spearman
   };

   ~Similarity();

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
   virtual QByteArray buildMPIBlock() override final;
   virtual bool readMPIBlock(const QByteArray& block) override final;
   virtual QByteArray processMPIBlock(const QByteArray& block) override final;

private:
   static const char* None;
   static const char* GMM;
   static const char* KMeans;
   static const char* Pearson;
   static const char* Spearman;
   static const char* BIC;
   static const char* ICL;

   int fetchPair(Pairwise::Index index, QVector<Pairwise::Vector2>& X, QVector<qint8>& labels);
   void savePair(Pairwise::Index index, qint8 K, const qint8 *labels, int N, const float *correlations);

   struct Block
   {
      enum
      {
         Start
         ,Load
         ,Execute1
         ,Execute2
         ,Execute3
         ,Read
         ,Done
      };

      Block(EOpenCLDevice& device, int N, int N_pow2, int K, int kernelSize)
      {
         in_index = device.makeBuffer<cl_int2>(1 * kernelSize).release();
         work_X = device.makeBuffer<Pairwise::Vector2>(N * kernelSize).release();
         work_N = device.makeBuffer<cl_int>(1 * kernelSize).release();
         work_labels = device.makeBuffer<cl_char>(N * kernelSize).release();
         work_components = device.makeBuffer<Pairwise::GMM::Component>(K * kernelSize).release();
         work_MP = device.makeBuffer<Pairwise::Vector2>(K * kernelSize).release();
         work_counts = device.makeBuffer<cl_int>(K * kernelSize).release();
         work_logpi = device.makeBuffer<cl_float>(K * kernelSize).release();
         work_loggamma = device.makeBuffer<cl_float>(N * K * kernelSize).release();
         work_logGamma = device.makeBuffer<cl_float>(K * kernelSize).release();
         out_K = device.makeBuffer<cl_char>(1 * kernelSize).release();
         out_labels = device.makeBuffer<cl_char>(N * kernelSize).release();

         work_x = device.makeBuffer<cl_float>(N_pow2 * kernelSize).release();
         work_y = device.makeBuffer<cl_float>(N_pow2 * kernelSize).release();
         work_rank = device.makeBuffer<cl_int>(N_pow2 * kernelSize).release();
         out_correlations = device.makeBuffer<cl_float>(K * kernelSize).release();

         if ( !device )
         {
            E_MAKE_EXCEPTION(e);
            throw e;
         }
      }

      ~Block()
      {
         delete in_index;
         delete work_X;
         delete work_N;
         delete work_labels;
         delete work_components;
         delete work_MP;
         delete work_counts;
         delete work_logpi;
         delete work_loggamma;
         delete work_logGamma;
         delete out_K;
         delete out_labels;

         delete work_x;
         delete work_y;
         delete work_rank;
         delete out_correlations;
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
      Pairwise::Index index;
      QVector<EOpenCLEvent> events;

      // clustering buffers
      EOpenCLBuffer<cl_int2>* in_index;
      EOpenCLBuffer<Pairwise::Vector2>* work_X;
      EOpenCLBuffer<cl_int>* work_N;
      EOpenCLBuffer<cl_char>* work_labels;
      EOpenCLBuffer<Pairwise::GMM::Component>* work_components;
      EOpenCLBuffer<Pairwise::Vector2>* work_MP;
      EOpenCLBuffer<cl_int>* work_counts;
      EOpenCLBuffer<cl_float>* work_logpi;
      EOpenCLBuffer<cl_float>* work_loggamma;
      EOpenCLBuffer<cl_float>* work_logGamma;
      EOpenCLBuffer<cl_char>* out_K;
      EOpenCLBuffer<cl_char>* out_labels;

      // correlation buffers
      EOpenCLBuffer<cl_float>* work_x;
      EOpenCLBuffer<cl_float>* work_y;
      EOpenCLBuffer<cl_int>* work_rank;
      EOpenCLBuffer<cl_float>* out_correlations;
   };

   void initializeOpenCL();
   void runStartBlock(Block& block);
   void runLoadBlock(Block& block);
   void runExecute1Block(Block& block);
   void runExecute2Block(Block& block);
   void runExecute3Block(Block& block);
   void runReadBlock(Block& block);

   ExpressionMatrix* _input {nullptr};
   CCMatrix* _ccm {nullptr};
   CorrelationMatrix* _cmx {nullptr};
   ClusteringMethod _clusMethod {ClusteringMethod::GMM};
   CorrelationMethod _corrMethod {CorrelationMethod::Pearson};
   Pairwise::Clustering* _clusModel {nullptr};
   Pairwise::Correlation* _corrModel {nullptr};
   int _minSamples {30};
   float _minExpression {-INFINITY};
   qint8 _minClusters {1};
   qint8 _maxClusters {5};
   Pairwise::Criterion _criterion {Pairwise::Criterion::ICL};
   bool _removePreOutliers {false};
   bool _removePostOutliers {false};
   float _minCorrelation {0.5};
   float _maxCorrelation {1.0};
   int _blockSize {4};
   int _kernelSize {4096};

   Block** _blocks {nullptr};
   EOpenCLProgram* _program {nullptr};
   EOpenCLKernel* _kernel1 {nullptr};
   EOpenCLKernel* _kernel2 {nullptr};
   EOpenCLKernel* _kernel3 {nullptr};
   EOpenCLBuffer<cl_float>* _expressions {nullptr};

   QDataStream *_mpiOut {nullptr};

   Pairwise::Index _index;
   Pairwise::Index _nextIndex;
   qint64 _totalSteps;
   qint64 _stepsStarted {0};
   qint64 _stepsComplete {0};
   int _lastPercent {0};
};



#endif
