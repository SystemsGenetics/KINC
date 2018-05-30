#include "similarity_opencl_kmeans.h"



using namespace std;






Similarity::OpenCL::KMeans::KMeans(::OpenCL::Program* program, QObject* parent):
   ::OpenCL::Kernel(program, "KMeans_compute", parent)
{
}






::OpenCL::Event Similarity::OpenCL::KMeans::execute(
   ::OpenCL::CommandQueue* queue,
   int kernelSize,
   ::OpenCL::Buffer<cl_float>* expressions,
   cl_int sampleSize,
   cl_int minSamples,
   cl_char minClusters,
   cl_char maxClusters,
   cl_int removePreOutliers,
   cl_int removePostOutliers,
   ::OpenCL::Buffer<Pairwise::Vector2>* work_X,
   ::OpenCL::Buffer<cl_int>* work_N,
   ::OpenCL::Buffer<cl_float>* work_outlier,
   ::OpenCL::Buffer<cl_char>* work_labels,
   ::OpenCL::Buffer<Pairwise::Vector2>* work_means,
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
   setArgument(RemovePreOutliers, removePreOutliers);
   setArgument(RemovePostOutliers, removePostOutliers);
   setBuffer(WorkX, work_X);
   setBuffer(WorkN, work_N);
   setBuffer(WorkOutlier, work_outlier);
   setBuffer(WorkLabels, work_labels);
   setBuffer(WorkMeans, work_means);
   setBuffer(OutK, out_K);
   setBuffer(OutLabels, out_labels);

   // set kernel sizes
   int workgroupSize {min(kernelSize, maxWorkGroupSize(queue->device())) / workGroupMultiple(queue->device()) * workGroupMultiple(queue->device())};
   setSizes(0, kernelSize, workgroupSize);

   // execute kernel
   return ::OpenCL::Kernel::execute(queue);
}
