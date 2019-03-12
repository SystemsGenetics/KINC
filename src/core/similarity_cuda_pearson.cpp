#include "similarity_cuda_pearson.h"



using namespace std;






/*!
 * Construct a new Pearson kernel object with the given CUDA program.
 *
 * @param program
 */
Similarity::CUDA::Pearson::Pearson(::CUDA::Program* program):
   ::CUDA::Kernel(program, "Pearson_compute")
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
 * @param out_correlations
 */
::CUDA::Event Similarity::CUDA::Pearson::execute(
   const ::CUDA::Stream& stream,
   int globalWorkSize,
   int localWorkSize,
   ::CUDA::Buffer<float2>* in_data,
   char clusterSize,
   ::CUDA::Buffer<qint8>* in_labels,
   int sampleSize,
   int minSamples,
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
      out_correlations);

   // set kernel arguments
   setArgument(GlobalWorkSize, globalWorkSize);
   setBuffer(InData, in_data);
   setArgument(ClusterSize, clusterSize);
   setBuffer(InLabels, in_labels);
   setArgument(SampleSize, sampleSize);
   setArgument(MinSamples, minSamples);
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
