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
 * @param numPairs
 * @param expressions
 * @param sampleSize
 * @param in_index
 * @param in_N
 * @param in_labels
 * @param in_K
 * @param marker
 * @param work_x
 * @param work_y
 */
::CUDA::Event Similarity::CUDA::Outlier::execute(
   const ::CUDA::Stream& stream,
   int globalWorkSize,
   int localWorkSize,
   int numPairs,
   ::CUDA::Buffer<float>* expressions,
   int sampleSize,
   ::CUDA::Buffer<int2>* in_index,
   ::CUDA::Buffer<int>* in_N,
   ::CUDA::Buffer<qint8>* in_labels,
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
      numPairs,
      expressions,
      sampleSize,
      in_index,
      in_N,
      in_labels,
      in_K,
      marker,
      work_x,
      work_y);

   // set kernel arguments
   setArgument(NumPairs, numPairs);
   setBuffer(Expressions, expressions);
   setArgument(SampleSize, sampleSize);
   setBuffer(InIndex, in_index);
   setBuffer(InN, in_N);
   setBuffer(InLabels, in_labels);
   setBuffer(InK, in_K);
   setArgument(Marker, marker);
   setBuffer(WorkX, work_x);
   setBuffer(WorkY, work_y);

   // set work sizes
   setSizes(globalWorkSize / localWorkSize, localWorkSize);

   // execute kernel
   return ::CUDA::Kernel::execute(stream);
}
