#include "similarity_cuda_spearman.h"



using namespace std;






/*!
 * Construct a new Spearman kernel object with the given CUDA program.
 *
 * @param program
 */
Similarity::CUDA::Spearman::Spearman(::CUDA::Program* program):
   ::CUDA::Kernel(program, "Spearman_compute")
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
 * @param in_data
 * @param clusterSize
 * @param in_labels
 * @param sampleSize
 * @param minSamples
 * @param work_xy
 * @param work_rank
 * @param out_correlations
 */
::CUDA::Event Similarity::CUDA::Spearman::execute(
   const ::CUDA::Stream& stream,
   int globalWorkSize,
   int localWorkSize,
   ::CUDA::Buffer<float2>* in_data,
   char clusterSize,
   ::CUDA::Buffer<qint8>* in_labels,
   int sampleSize,
   int minSamples,
   ::CUDA::Buffer<float>* work_xy,
   ::CUDA::Buffer<int>* work_rank,
   ::CUDA::Buffer<float>* out_correlations
)
{
   EDEBUG_FUNC(this,
      &stream,
      globalWorkSize,
      localWorkSize,
      in_data,
      clusterSize,
      in_labels,
      sampleSize,
      minSamples,
      work_xy,
      work_rank,
      out_correlations);

   // set kernel arguments
   setArgument(GlobalWorkSize, globalWorkSize);
   setBuffer(InData, in_data);
   setArgument(ClusterSize, clusterSize);
   setBuffer(InLabels, in_labels);
   setArgument(SampleSize, sampleSize);
   setArgument(MinSamples, minSamples);
   setBuffer(WorkXY, work_xy);
   setBuffer(WorkRank, work_rank);
   setBuffer(OutCorrelations, out_correlations);

   // set work sizes
   if ( localWorkSize == 0 )
   {
      localWorkSize = min(globalWorkSize, getAttribute(CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK));
   }

   int numWorkgroups = (globalWorkSize + localWorkSize - 1) / localWorkSize;

   setSizes(numWorkgroups * localWorkSize, localWorkSize);

   // execute kernel
   return ::CUDA::Kernel::execute(stream);
}
