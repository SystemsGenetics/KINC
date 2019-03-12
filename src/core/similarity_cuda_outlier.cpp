#include "similarity_cuda_outlier.h"



using namespace std;






/*!
 * Construct a new Outlier kernel object with the given OpenCL program.
 *
 * @param program
 */
Similarity::CUDA::Outlier::Outlier(::CUDA::Program* program):
   ::CUDA::Kernel(program, "removeOutliers")
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
 * @param in_N
 * @param in_labels
 * @param sampleSize
 * @param in_K
 * @param marker
 * @param work_x
 * @param work_y
 */
::CUDA::Event Similarity::CUDA::Outlier::execute(
   const ::CUDA::Stream& stream,
   int globalWorkSize,
   int localWorkSize,
   ::CUDA::Buffer<float2>* in_data,
   ::CUDA::Buffer<int>* in_N,
   ::CUDA::Buffer<qint8>* in_labels,
   int sampleSize,
   ::CUDA::Buffer<qint8>* in_K,
   qint8 marker,
   ::CUDA::Buffer<float>* work_x,
   ::CUDA::Buffer<float>* work_y
)
{
   EDEBUG_FUNC(this,
      &stream,
      globalWorkSize,
      localWorkSize,
      in_data,
      in_N,
      in_labels,
      sampleSize,
      in_K,
      marker,
      work_x,
      work_y);

   // set kernel arguments
   setArgument(GlobalWorkSize, globalWorkSize);
   setBuffer(InData, in_data);
   setBuffer(InN, in_N);
   setBuffer(InLabels, in_labels);
   setArgument(SampleSize, sampleSize);
   setBuffer(InK, in_K);
   setArgument(Marker, marker);
   setBuffer(WorkX, work_x);
   setBuffer(WorkY, work_y);

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
