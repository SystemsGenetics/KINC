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
 * @param work_gmm_pi
 * @param work_gmm_mu
 * @param work_gmm_sigma
 * @param work_gmm_sigmaInv
 * @param work_gmm_normalizer
 * @param work_gmm_MP
 * @param work_gmm_counts
 * @param work_gmm_logpi
 * @param work_gmm_gamma
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
   ::CUDA::Buffer<float>* work_gmm_pi,
   ::CUDA::Buffer<float2>* work_gmm_mu,
   ::CUDA::Buffer<float4>* work_gmm_sigma,
   ::CUDA::Buffer<float4>* work_gmm_sigmaInv,
   ::CUDA::Buffer<float>* work_gmm_normalizer,
   ::CUDA::Buffer<float2>* work_gmm_MP,
   ::CUDA::Buffer<int>* work_gmm_counts,
   ::CUDA::Buffer<float>* work_gmm_logpi,
   ::CUDA::Buffer<float>* work_gmm_gamma,
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
      work_gmm_pi,
      work_gmm_mu,
      work_gmm_sigma,
      work_gmm_sigmaInv,
      work_gmm_normalizer,
      work_gmm_MP,
      work_gmm_counts,
      work_gmm_logpi,
      work_gmm_gamma,
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
   setBuffer(WorkGmmPi, work_gmm_pi);
   setBuffer(WorkGmmMu, work_gmm_mu);
   setBuffer(WorkGmmSigma, work_gmm_sigma);
   setBuffer(WorkGmmSigmaInv, work_gmm_sigmaInv);
   setBuffer(WorkGmmNormalizer, work_gmm_normalizer);
   setBuffer(WorkGmmMP, work_gmm_MP);
   setBuffer(WorkGmmCounts, work_gmm_counts);
   setBuffer(WorkGmmLogPi, work_gmm_logpi);
   setBuffer(WorkGmmGamma, work_gmm_gamma);
   setBuffer(OutK, out_K);
   setBuffer(OutLabels, out_labels);

   // set work sizes
   setSizes(globalWorkSize / localWorkSize, localWorkSize);

   // execute kernel
   return ::CUDA::Kernel::execute(stream);
}
