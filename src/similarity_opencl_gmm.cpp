#include "similarity_opencl_gmm.h"



using namespace std;






Similarity::OpenCL::GMM::GMM(::OpenCL::Program* program, QObject* parent):
   ::OpenCL::Kernel(program, "GMM_compute", parent)
{
}






::OpenCL::Event Similarity::OpenCL::GMM::execute(
   ::OpenCL::CommandQueue* queue,
   int kernelSize,
   ::OpenCL::Buffer<cl_float>* expressions,
   cl_int sampleSize,
   cl_int minSamples,
   cl_char minClusters,
   cl_char maxClusters,
   Pairwise::Criterion criterion,
   cl_int removePreOutliers,
   cl_int removePostOutliers,
   ::OpenCL::Buffer<Pairwise::Vector2>* work_X,
   ::OpenCL::Buffer<cl_int>* work_N,
   ::OpenCL::Buffer<cl_char>* work_labels,
   ::OpenCL::Buffer<Pairwise::GMM::Component>* work_components,
   ::OpenCL::Buffer<Pairwise::Vector2>* work_MP,
   ::OpenCL::Buffer<cl_int>* work_counts,
   ::OpenCL::Buffer<cl_float>* work_logpi,
   ::OpenCL::Buffer<cl_float>* work_loggamma,
   ::OpenCL::Buffer<cl_float>* work_logGamma,
   ::OpenCL::Buffer<cl_char>* out_K,
   ::OpenCL::Buffer<cl_char>* out_labels
)
{
   // acquire lock for this kernel
   Locker locker {lock()};

   // set kernel arguments
   setBuffer(Expressions, expressions);
   setArgument(SampleSize, sampleSize);
   setArgument(MinSamples, minSamples);
   setArgument(MinClusters, minClusters);
   setArgument(MaxClusters, maxClusters);
   setArgument(Criterion, criterion);
   setArgument(RemovePreOutliers, removePreOutliers);
   setArgument(RemovePostOutliers, removePostOutliers);
   setBuffer(WorkX, work_X);
   setBuffer(WorkN, work_N);
   setBuffer(WorkLabels, work_labels);
   setBuffer(WorkComponents, work_components);
   setBuffer(WorkMP, work_MP);
   setBuffer(WorkCounts, work_counts);
   setBuffer(WorkLogPi, work_logpi);
   setBuffer(WorkLoggamma, work_loggamma);
   setBuffer(WorkLogGamma, work_logGamma);
   setBuffer(OutK, out_K);
   setBuffer(OutLabels, out_labels);

   // set kernel sizes
   setSizes(0, kernelSize, min(kernelSize, maxWorkGroupSize(queue->device())));

   // execute kernel
   return ::OpenCL::Kernel::execute(queue);
}
