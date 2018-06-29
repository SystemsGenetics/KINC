#ifndef SIMILARITY_OPENCL_SPEARMAN_H
#define SIMILARITY_OPENCL_SPEARMAN_H
#include "similarity_opencl.h"



class Similarity::OpenCL::Spearman : public ::OpenCL::Kernel
{
   Q_OBJECT
public:
   enum Argument
   {
      InData
      ,ClusterSize
      ,InLabels
      ,SampleSize
      ,MinSamples
      ,WorkX
      ,WorkY
      ,WorkRank
      ,OutCorrelations
   };
   explicit Spearman(::OpenCL::Program* program, QObject* parent = nullptr);
   ::OpenCL::Event execute(
      ::OpenCL::CommandQueue* queue,
      int kernelSize,
      ::OpenCL::Buffer<Pairwise::Vector2>* in_data,
      cl_char clusterSize,
      ::OpenCL::Buffer<cl_char>* in_labels,
      cl_int sampleSize,
      cl_int minSamples,
      ::OpenCL::Buffer<cl_float>* work_x,
      ::OpenCL::Buffer<cl_float>* work_y,
      ::OpenCL::Buffer<cl_int>* work_rank,
      ::OpenCL::Buffer<cl_float>* out_correlations
   );
};



#endif
