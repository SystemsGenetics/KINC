#include "similarity_cuda_fetchpair.h"



using namespace std;






/*!
 * Construct a new fetch-pair kernel object with the given OpenCL program.
 *
 * @param program
 */
Similarity::CUDA::FetchPair::FetchPair(::CUDA::Program* program):
   ::CUDA::Kernel(program, "fetchPair")
{
   EDEBUG_FUNC(this,program);
}






/*!
 * Execute this kernel object's OpenCL kernel using the given OpenCL command
 * queue and kernel arguments, returning the OpenCL event associated with the
 * kernel execution.
 *
 * @param queue
 * @param globalWorkSize
 * @param localWorkSize
 * @param expressions
 * @param sampleSize
 * @param in_index
 * @param minExpression
 * @param out_X
 * @param out_N
 * @param out_labels
 */
::CUDA::Event Similarity::CUDA::FetchPair::execute(
   const ::CUDA::Stream& stream,
   int globalWorkSize,
   int localWorkSize,
   ::CUDA::Buffer<float>* expressions,
   int sampleSize,
   ::CUDA::Buffer<int2>* in_index,
   int minExpression,
   ::CUDA::Buffer<float2>* out_X,
   ::CUDA::Buffer<int>* out_N,
   ::CUDA::Buffer<qint8>* out_labels
)
{
   EDEBUG_FUNC(this,
      stream,
      globalWorkSize,
      localWorkSize,
      expressions,
      sampleSize,
      in_index,
      minExpression,
      out_X,
      out_N,
      out_labels);

   // set kernel arguments
   setArgument(GlobalWorkSize, globalWorkSize);
   setBuffer(Expressions, expressions);
   setArgument(SampleSize, sampleSize);
   setBuffer(InIndex, in_index);
   setArgument(MinExpression, minExpression);
   setBuffer(OutX, out_X);
   setBuffer(OutN, out_N);
   setBuffer(OutLabels, out_labels);

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
