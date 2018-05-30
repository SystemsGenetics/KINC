#include "similarity_opencl_pearson.h"



using namespace std;






Similarity::OpenCL::Pearson::Pearson(::OpenCL::Program* program, QObject* parent):
   ::OpenCL::Kernel(program, "Pearson_compute", parent)
{
}






::OpenCL::Event Similarity::OpenCL::Pearson::execute(
   ::OpenCL::CommandQueue* queue,
   int kernelSize,
   ::OpenCL::Buffer<Pairwise::Vector2>* in_data,
   cl_char clusterSize,
   ::OpenCL::Buffer<cl_char>* in_labels,
   cl_int sampleSize,
   cl_int minSamples,
   ::OpenCL::Buffer<cl_float>* out_correlations
)
{
   // acquire lock for this kernel
   Locker locker {lock()};

   // set kernel arguments
   setBuffer(InData, in_data);
   setArgument(ClusterSize, clusterSize);
   setBuffer(InLabels, in_labels);
   setArgument(SampleSize, sampleSize);
   setArgument(MinSamples, minSamples);
   setBuffer(OutCorrelations, out_correlations);

   // set kernel sizes
   setSizes(0, kernelSize, min(kernelSize, maxWorkGroupSize(queue->device())));

   // execute kernel
   return ::OpenCL::Kernel::execute(queue);
}
