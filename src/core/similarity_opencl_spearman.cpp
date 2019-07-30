#include "similarity_opencl_spearman.h"



using namespace std;






/*!
 * Construct a new Spearman kernel object with the given OpenCL program and
 * qt parent.
 *
 * @param program
 * @param parent
 */
Similarity::OpenCL::Spearman::Spearman(::OpenCL::Program* program, QObject* parent):
   ::OpenCL::Kernel(program, "Spearman_compute", parent)
{
   EDEBUG_FUNC(this,parent);
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
 * @param clusterSize
 * @param in_labels
 * @param minSamples
 * @param work_x
 * @param work_y
 * @param out_correlations
 */
::OpenCL::Event Similarity::OpenCL::Spearman::execute(
   ::OpenCL::CommandQueue* queue,
   int globalWorkSize,
   int localWorkSize,
   int numPairs,
   ::OpenCL::Buffer<cl_float>* expressions,
   cl_int sampleSize,
   ::OpenCL::Buffer<cl_int2>* in_index,
   cl_char clusterSize,
   ::OpenCL::Buffer<cl_char>* in_labels,
   cl_int minSamples,
   ::OpenCL::Buffer<cl_float>* work_x,
   ::OpenCL::Buffer<cl_float>* work_y,
   ::OpenCL::Buffer<cl_float>* out_correlations
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
      clusterSize,
      in_labels,
      minSamples,
      work_x,
      work_y,
      out_correlations);

   // acquire lock for this kernel
   Locker locker {lock()};

   // set kernel arguments
   setArgument(NumPairs, numPairs);
   setBuffer(Expressions, expressions);
   setArgument(SampleSize, sampleSize);
   setBuffer(InIndex, in_index);
   setArgument(ClusterSize, clusterSize);
   setBuffer(InLabels, in_labels);
   setArgument(MinSamples, minSamples);
   setBuffer(WorkX, work_x);
   setBuffer(WorkY, work_y);
   setBuffer(OutCorrelations, out_correlations);

   // set work sizes
   setSizes(0, globalWorkSize / localWorkSize, localWorkSize);

   // execute kernel
   return ::OpenCL::Kernel::execute(queue);
}
