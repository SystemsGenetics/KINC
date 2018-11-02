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
 * @param sampleSize
 * @param minSamples
 * @param minClusters
 * @param maxClusters
 * @param criterion
 * @param removePreOutliers
 * @param removePostOutliers
 * @param work_X
 * @param work_N
 * @param work_x
 * @param work_y
 * @param work_labels
 * @param work_components
 * @param work_MP
 * @param work_counts
 * @param work_logpi
 * @param work_loggamma
 * @param work_logGamma
 * @param out_K
 * @param out_labels
 */
::OpenCL::Event Similarity::OpenCL::GMM::execute(
   ::OpenCL::CommandQueue* queue,
   int globalWorkSize,
   int localWorkSize,
   cl_int sampleSize,
   cl_int minSamples,
   cl_char minClusters,
   cl_char maxClusters,
   cl_int criterion,
   cl_int removePreOutliers,
   cl_int removePostOutliers,
   ::OpenCL::Buffer<cl_float2>* work_X,
   ::OpenCL::Buffer<cl_int>* work_N,
   ::OpenCL::Buffer<cl_float>* work_x,
   ::OpenCL::Buffer<cl_float>* work_y,
   ::OpenCL::Buffer<cl_char>* work_labels,
   ::OpenCL::Buffer<cl_component>* work_components,
   ::OpenCL::Buffer<cl_float2>* work_MP,
   ::OpenCL::Buffer<cl_int>* work_counts,
   ::OpenCL::Buffer<cl_float>* work_logpi,
   ::OpenCL::Buffer<cl_float>* work_loggamma,
   ::OpenCL::Buffer<cl_float>* work_logGamma,
   ::OpenCL::Buffer<cl_char>* out_K,
   ::OpenCL::Buffer<cl_char>* out_labels
)
{
   EDEBUG_FUNC(this,
      queue,
      globalWorkSize,
      localWorkSize,
      sampleSize,
      minSamples,
      minClusters,
      maxClusters,
      &criterion,
      removePreOutliers,
      removePostOutliers,
      work_X,
      work_N,
      work_x,
      work_y,
      work_labels,
      work_components,
      work_MP,
      work_counts,
      work_logpi,
      work_loggamma,
      work_logGamma,
      out_K,
      out_labels);

   // acquire lock for this kernel
   Locker locker {lock()};

   // set kernel arguments
   setArgument(GlobalWorkSize, globalWorkSize);
   setArgument(SampleSize, sampleSize);
   setArgument(MinSamples, minSamples);
   setArgument(MinClusters, minClusters);
   setArgument(MaxClusters, maxClusters);
   setArgument(Criterion, criterion);
   setArgument(RemovePreOutliers, removePreOutliers);
   setArgument(RemovePostOutliers, removePostOutliers);
   setBuffer(WorkX, work_X);
   setBuffer(WorkN, work_N);
   setBuffer(WorkXSorted, work_x);
   setBuffer(WorkYSorted, work_y);
   setBuffer(WorkLabels, work_labels);
   setBuffer(WorkComponents, work_components);
   setBuffer(WorkMP, work_MP);
   setBuffer(WorkCounts, work_counts);
   setBuffer(WorkLogPi, work_logpi);
   setBuffer(WorkLoggamma, work_loggamma);
   setBuffer(WorkLogGamma, work_logGamma);
   setBuffer(OutK, out_K);
   setBuffer(OutLabels, out_labels);

   // set work sizes
   if ( localWorkSize == 0 )
   {
      localWorkSize = min(globalWorkSize, maxWorkGroupSize(queue->device()));
   }

   int numWorkgroups = (globalWorkSize + localWorkSize - 1) / localWorkSize;

   setSizes(0, numWorkgroups * localWorkSize, localWorkSize);

   // execute kernel
   return ::OpenCL::Kernel::execute(queue);
}
