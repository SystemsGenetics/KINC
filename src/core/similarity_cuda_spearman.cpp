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
 * @param expressions
 * @param sampleSize
 * @param in_index
 * @param clusterSize
 * @param in_labels
 * @param minSamples
 * @param work_xy
 * @param out_correlations
 */
::CUDA::Event Similarity::CUDA::Spearman::execute(
   const ::CUDA::Stream& stream,
   int globalWorkSize,
   int localWorkSize,
   ::CUDA::Buffer<float>* expressions,
   int sampleSize,
   ::CUDA::Buffer<int2>* in_index,
   char clusterSize,
   ::CUDA::Buffer<qint8>* in_labels,
   int minSamples,
   ::CUDA::Buffer<float>* work_xy,
   ::CUDA::Buffer<float>* out_correlations
)
{
   EDEBUG_FUNC(this,
      &stream,
      globalWorkSize,
      localWorkSize,
      expressions,
      sampleSize,
      in_index,
      clusterSize,
      in_labels,
      minSamples,
      work_xy,
      out_correlations);

   // set kernel arguments
   setArgument(GlobalWorkSize, globalWorkSize);
   setBuffer(Expressions, expressions);
   setArgument(SampleSize, sampleSize);
   setBuffer(InIndex, in_index);
   setArgument(ClusterSize, clusterSize);
   setBuffer(InLabels, in_labels);
   setArgument(MinSamples, minSamples);
   setBuffer(WorkXY, work_xy);
   setBuffer(OutCorrelations, out_correlations);

   // set work sizes
   int numWorkgroups = (globalWorkSize + localWorkSize - 1) / localWorkSize;

   setSizes(numWorkgroups, localWorkSize);

   // execute kernel
   return ::CUDA::Kernel::execute(stream);
}
