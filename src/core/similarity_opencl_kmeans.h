#ifndef SIMILARITY_OPENCL_KMEANS_H
#define SIMILARITY_OPENCL_KMEANS_H
#include "similarity_opencl.h"



class Similarity::OpenCL::KMeans : public ::OpenCL::Kernel
{
   Q_OBJECT
public:
   enum Argument
   {
      Expressions
      ,SampleSize
      ,MinSamples
      ,MinClusters
      ,MaxClusters
      ,RemovePreOutliers
      ,RemovePostOutliers
      ,WorkX
      ,WorkN
      ,WorkOutlier
      ,WorkLabels
      ,WorkMeans
      ,OutK
      ,OutLabels
   };
   explicit KMeans(::OpenCL::Program* program, QObject* parent = nullptr);
   ::OpenCL::Event execute(
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
   );
};



#endif
