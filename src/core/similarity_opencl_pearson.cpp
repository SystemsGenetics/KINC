#include "similarity_opencl_pearson.h"



using namespace std;






/*!
 * Construct a new Pearson kernel object with the given OpenCL program and
 * qt parent.
 *
 * @param program
 * @param parent
 */
Similarity::OpenCL::Pearson::Pearson(::OpenCL::Program* program, QObject* parent):
   ::OpenCL::Kernel(program, "Pearson_compute", parent)
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
 * @param in_data
 * @param clusterSize
 * @param in_labels
 * @param sampleSize
 * @param minSamples
 * @param out_correlations
 */
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
   EDEBUG_FUNC(this,
      queue,
      kernelSize,
      in_data,
      clusterSize,
      in_labels,
      sampleSize,
      minSamples,
      out_correlations);

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
   int localSize {min(kernelSize, maxWorkGroupSize(queue->device()))};
   int globalSize {kernelSize - kernelSize % localSize};

   setSizes(0, globalSize, localSize);

   // execute kernel
   return ::OpenCL::Kernel::execute(queue);
}
