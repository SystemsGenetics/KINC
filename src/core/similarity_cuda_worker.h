#ifndef SIMILARITY_CUDA_WORKER_H
#define SIMILARITY_CUDA_WORKER_H
#include "similarity_cuda.h"
#include "similarity_cuda_fetchpair.h"
#include "similarity_cuda_gmm.h"
#include "similarity_cuda_outlier.h"
#include "similarity_cuda_pearson.h"
#include "similarity_cuda_spearman.h"



/*!
 * This class implements the CUDA worker of the similarity analytic.
 */
class Similarity::CUDA::Worker : public EAbstractAnalyticCUDAWorker
{
   Q_OBJECT
public:
   explicit Worker(Similarity* base, Similarity::CUDA* baseCuda, ::CUDA::Program* program);
   virtual std::unique_ptr<EAbstractAnalyticBlock> execute(const EAbstractAnalyticBlock* block) override final;
private:
   /*!
    * Pointer to the base analytic.
    */
   Similarity* _base;
   /*!
    * Pointer to the base CUDA object.
    */
   Similarity::CUDA* _baseCuda;
   /*!
    * Pointer to this worker's unique and private stream.
    */
   ::CUDA::Stream _stream;
   /*!
    * Structure of this worker's kernels.
    */
   struct
   {
      CUDA::FetchPair fetchPair;
      CUDA::GMM gmm;
      CUDA::Outlier outlier;
      CUDA::Pearson pearson;
      CUDA::Spearman spearman;
   } _kernels;
   /*!
    * Structure of this worker's buffers.
    */
   struct
   {
      ::CUDA::Buffer<int2> in_index;
      ::CUDA::Buffer<int> work_N;
      ::CUDA::Buffer<float> work_xy;
      ::CUDA::Buffer<qint8> work_labels;
      ::CUDA::Buffer<cu_component> work_components;
      ::CUDA::Buffer<float2> work_MP;
      ::CUDA::Buffer<int> work_counts;
      ::CUDA::Buffer<float> work_logpi;
      ::CUDA::Buffer<float> work_gamma;
      ::CUDA::Buffer<qint8> out_K;
      ::CUDA::Buffer<qint8> out_labels;
      ::CUDA::Buffer<float> out_correlations;
   } _buffers;
};



#endif
