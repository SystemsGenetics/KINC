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
::OpenCL::Event Similarity::OpenCL::FetchPair::execute(
   ::OpenCL::CommandQueue* queue,
   int globalWorkSize,
   int localWorkSize,
   int numPairs,
   ::OpenCL::Buffer<cl_float>* expressions,
   cl_int sampleSize,
   ::OpenCL::Buffer<cl_int2>* in_index,
   cl_float minExpression,
   ::OpenCL::Buffer<cl_int>* out_N,
   ::OpenCL::Buffer<cl_char>* out_labels
)
{
   EDEBUG_FUNC(this,
      queue,
      globalWorkSize,
      localWorkSize,
      numPairs,
      expressions,
      sampleSize,
      in_index,
      minExpression,
      out_N,
      out_labels);

   // acquire lock for this kernel
   Locker locker {lock()};

   // set kernel arguments
   setArgument(NumPairs, numPairs);
   setBuffer(Expressions, expressions);
   setArgument(SampleSize, sampleSize);
   setBuffer(InIndex, in_index);
   setArgument(MinExpression, minExpression);
   setBuffer(OutN, out_N);
   setBuffer(OutLabels, out_labels);

   // set work sizes
   setSizes(0, globalWorkSize / localWorkSize, localWorkSize);

   // execute kernel
   return ::OpenCL::Kernel::execute(queue);
}
