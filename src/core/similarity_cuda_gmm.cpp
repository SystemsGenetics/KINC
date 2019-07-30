#include "similarity_cuda_gmm.h"



using namespace std;






/*!
 * Construct a new GMM kernel object with the given CUDA program.
 *
 * @param program
 */
Similarity::CUDA::GMM::GMM(::CUDA::Program* program):
   ::CUDA::Kernel(program, "GMM_compute")
{
   EDEBUG_FUNC(this,program);
}






/*!
 * Execute this kernel object's CUDA kernel using the given CUDA stream
 * and kernel arguments, returning the CUDA event associated with the
 * kernel execution.
 *
 * @param stream
 * @param globalWorkSize
 * @param localWorkSize
 * @param numPairs
 * @param expressions
 * @param sampleSize
 * @param in_index
 * @param minSamples
 * @param minClusters
 * @param maxClusters
 * @param criterion
 * @param work_X
 * @param work_N
 * @param work_labels
 * @param work_components
 * @param work_MP
 * @param work_counts
 * @param work_logpi
 * @param work_gamma
 * @param out_K
 * @param out_labels
 */
::CUDA::Event Similarity::CUDA::GMM::execute(
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
)
{
   EDEBUG_FUNC(this,
      &stream,
      globalWorkSize,
      localWorkSize,
      numPairs,
      expressions,
      sampleSize,
      in_index,
      minSamples,
      minClusters,
      maxClusters,
      criterion,
      work_X,
      work_N,
      work_labels,
      work_components,
      work_MP,
      work_counts,
      work_logpi,
      work_gamma,
      out_K,
      out_labels);

   // set kernel arguments
   setArgument(NumPairs, numPairs);
   setBuffer(Expressions, expressions);
   setArgument(SampleSize, sampleSize);
   setBuffer(InIndex, in_index);
   setArgument(MinSamples, minSamples);
   setArgument(MinClusters, minClusters);
   setArgument(MaxClusters, maxClusters);
   setArgument(Criterion, criterion);
   setBuffer(WorkX, work_X);
   setBuffer(WorkN, work_N);
   setBuffer(WorkLabels, work_labels);
   setBuffer(WorkComponents, work_components);
   setBuffer(WorkMP, work_MP);
   setBuffer(WorkCounts, work_counts);
   setBuffer(WorkLogPi, work_logpi);
   setBuffer(WorkGamma, work_gamma);
   setBuffer(OutK, out_K);
   setBuffer(OutLabels, out_labels);

   // set work sizes
   setSizes(globalWorkSize / localWorkSize, localWorkSize);

   // execute kernel
   return ::CUDA::Kernel::execute(stream);
}
