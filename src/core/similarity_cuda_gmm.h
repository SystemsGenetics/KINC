#ifndef SIMILARITY_CUDA_GMM_H
#define SIMILARITY_CUDA_GMM_H
#include "similarity_cuda.h"



typedef struct
{
   float pi;
   float2 mu;
   float4 sigma;
   float4 sigmaInv;
   float normalizer;
} cu_component;






/*!
 * This class implements the GMM kernel for the similarity analytic. This
 * kernel takes a list of pairwise data arrays and computes the number of
 * clusters and a list of cluster labels for each pair.
 */
class Similarity::CUDA::GMM : public ::CUDA::Kernel
{
public:
   /*!
    * Defines the arguments passed to the OpenCL kernel.
    */
   enum Argument
   {
      NumPairs
      ,Expressions
      ,SampleSize
      ,InIndex
      ,MinSamples
      ,MinClusters
      ,MaxClusters
      ,Criterion
      ,WorkX
      ,WorkN
      ,WorkLabels
      ,WorkComponents
      ,WorkMP
      ,WorkCounts
      ,WorkLogPi
      ,WorkGamma
      ,OutK
      ,OutLabels
   };
   explicit GMM(::CUDA::Program* program);
   ::CUDA::Event execute(
      const ::CUDA::Stream& stream,
      int globalWorkSize,
      int localWorkSize,
      int numPairs,
      ::CUDA::Buffer<float>* expressions,
      int sampleSize,
      ::CUDA::Buffer<int2>* in_index,
      int minSamples,
      char minClusters,
      char maxClusters,
      int criterion,
      ::CUDA::Buffer<float2>* work_X,
      ::CUDA::Buffer<int>* work_N,
      ::CUDA::Buffer<qint8>* work_labels,
      ::CUDA::Buffer<cu_component>* work_components,
      ::CUDA::Buffer<float2>* work_MP,
      ::CUDA::Buffer<int>* work_counts,
      ::CUDA::Buffer<float>* work_logpi,
      ::CUDA::Buffer<float>* work_gamma,
      ::CUDA::Buffer<qint8>* out_K,
      ::CUDA::Buffer<qint8>* out_labels
   );
};



#endif
