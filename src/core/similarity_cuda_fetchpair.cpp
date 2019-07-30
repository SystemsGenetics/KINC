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
 * @param numPairs
 * @param expressions
 * @param sampleSize
 * @param in_index
 * @param minExpression
 * @param out_N
 * @param out_labels
 */
::CUDA::Event Similarity::CUDA::FetchPair::execute(
   const ::CUDA::Stream& stream,
   int globalWorkSize,
   int localWorkSize,
   int numPairs,
   ::CUDA::Buffer<float>* expressions,
   int sampleSize,
   ::CUDA::Buffer<int2>* in_index,
   float minExpression,
   ::CUDA::Buffer<int>* out_N,
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
      minExpression,
      out_N,
      out_labels);

   // set kernel arguments
   setArgument(NumPairs, numPairs);
   setBuffer(Expressions, expressions);
   setArgument(SampleSize, sampleSize);
   setBuffer(InIndex, in_index);
   setArgument(MinExpression, minExpression);
   setBuffer(OutN, out_N);
   setBuffer(OutLabels, out_labels);

   // set work sizes
   setSizes(globalWorkSize / localWorkSize, localWorkSize);

   // execute kernel
   return ::CUDA::Kernel::execute(stream);
}
