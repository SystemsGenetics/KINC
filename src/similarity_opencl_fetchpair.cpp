#include "similarity_opencl_fetchpair.h"



using namespace std;






Similarity::OpenCL::FetchPair::FetchPair(::OpenCL::Program* program, QObject* parent):
   ::OpenCL::Kernel(program, "fetchPair", parent)
{
}






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
   int workgroupSize {min(kernelSize, maxWorkGroupSize(queue->device())) / workGroupMultiple(queue->device()) * workGroupMultiple(queue->device())};
   setSizes(0, kernelSize, workgroupSize);

   // execute kernel
   return ::OpenCL::Kernel::execute(queue);
}
