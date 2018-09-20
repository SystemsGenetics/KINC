#include "similarity_opencl_fetchpair.h"



using namespace std;






/*!
 * Construct a new fetch-pair kernel object with the given OpenCL program and
 * qt parent.
 *
 * @param program
 * @param parent
 */
Similarity::OpenCL::FetchPair::FetchPair(::OpenCL::Program* program, QObject* parent):
   ::OpenCL::Kernel(program, "fetchPair", parent)
{
   EDEBUG_FUNC(this,program,parent);
}






/*!
 * Execute this kernel object's OpenCL kernel using the given OpenCL command
 * queue and kernel arguments, returning the OpenCL event associated with the
 * kernel execution.
 *
 * @param queue
 * @param kernelSize
 * @param expressions
 * @param sampleSize
 * @param in_index
 * @param minExpression
 * @param out_X
 * @param out_N
 * @param out_labels
 */
::OpenCL::Event Similarity::OpenCL::FetchPair::execute(
   ::OpenCL::CommandQueue* queue,
   int kernelSize,
   ::OpenCL::Buffer<cl_float>* expressions,
   cl_int sampleSize,
   ::OpenCL::Buffer<cl_int2>* in_index,
   cl_int minExpression,
   ::OpenCL::Buffer<Pairwise::Vector2>* out_X,
   ::OpenCL::Buffer<cl_int>* out_N,
   ::OpenCL::Buffer<cl_char>* out_labels
)
{
   EDEBUG_FUNC(this,
      queue,
      kernelSize,
      expressions,
      sampleSize,
      in_index,
      minExpression,
      out_X,
      out_N,
      out_labels);

   // acquire lock for this kernel
   Locker locker {lock()};

   // set kernel arguments
   setBuffer(Expressions, expressions);
   setArgument(SampleSize, sampleSize);
   setBuffer(InIndex, in_index);
   setArgument(MinExpression, minExpression);
   setBuffer(OutX, out_X);
   setBuffer(OutN, out_N);
   setBuffer(OutLabels, out_labels);

   // set kernel sizes
   setSizes(0, kernelSize, min(kernelSize, maxWorkGroupSize(queue->device())));

   // execute kernel
   return ::OpenCL::Kernel::execute(queue);
}
