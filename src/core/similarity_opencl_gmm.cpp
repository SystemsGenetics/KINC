#include "similarity_opencl_gmm.h"



using namespace std;






/*!
 * Construct a new GMM kernel object with the given OpenCL program and qt parent.
 *
 * @param program
 * @param parent
 */
Similarity::OpenCL::GMM::GMM(::OpenCL::Program* program, QObject* parent):
   ::OpenCL::Kernel(program, "GMM_compute", parent)
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
 * @param minSamples
 * @param minClusters
 * @param maxClusters
 * @param criterion
 * @param work_X
 * @param work_N
 * @param work_labels
 * @param work_components
 * @param work_MP
 * @param work_counts
 * @param work_logpi
 * @param work_gamma
 * @param out_K
 * @param out_labels
 */
::OpenCL::Event Similarity::OpenCL::GMM::execute(
   ::OpenCL::CommandQueue* queue,
   int globalWorkSize,
   int localWorkSize,
   int numPairs,
   ::OpenCL::Buffer<cl_float>* expressions,
   cl_int sampleSize,
   ::OpenCL::Buffer<cl_int2>* in_index,
   cl_int minSamples,
   cl_char minClusters,
   cl_char maxClusters,
   cl_int criterion,
   ::OpenCL::Buffer<cl_float2>* work_X,
   ::OpenCL::Buffer<cl_int>* work_N,
   ::OpenCL::Buffer<cl_char>* work_labels,
   ::OpenCL::Buffer<cl_component>* work_components,
   ::OpenCL::Buffer<cl_float2>* work_MP,
   ::OpenCL::Buffer<cl_int>* work_counts,
   ::OpenCL::Buffer<cl_float>* work_logpi,
   ::OpenCL::Buffer<cl_float>* work_gamma,
   ::OpenCL::Buffer<cl_char>* out_K,
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
      minSamples,
      minClusters,
      maxClusters,
      &criterion,
      work_X,
      work_N,
      work_labels,
      work_components,
      work_MP,
      work_counts,
      work_logpi,
      work_gamma,
      out_K,
      out_labels);

   // acquire lock for this kernel
   Locker locker {lock()};

   // set kernel arguments
   setArgument(NumPairs, numPairs);
   setBuffer(Expressions, expressions);
   setArgument(SampleSize, sampleSize);
   setBuffer(InIndex, in_index);
   setArgument(MinSamples, minSamples);
   setArgument(MinClusters, minClusters);
   setArgument(MaxClusters, maxClusters);
   setArgument(Criterion, criterion);
   setBuffer(WorkX, work_X);
   setBuffer(WorkN, work_N);
   setBuffer(WorkLabels, work_labels);
   setBuffer(WorkComponents, work_components);
   setBuffer(WorkMP, work_MP);
   setBuffer(WorkCounts, work_counts);
   setBuffer(WorkLogPi, work_logpi);
   setBuffer(WorkGamma, work_gamma);
   setBuffer(OutK, out_K);
   setBuffer(OutLabels, out_labels);

   // set work sizes
   setSizes(0, globalWorkSize / localWorkSize, localWorkSize);

   // execute kernel
   return ::OpenCL::Kernel::execute(queue);
}
