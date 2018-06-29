#ifndef SIMILARITY_OPENCL_PEARSON_H
#define SIMILARITY_OPENCL_PEARSON_H
#include "similarity_opencl.h"



class Similarity::OpenCL::Pearson : public ::OpenCL::Kernel
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
      ,OutCorrelations
   };
   explicit Pearson(::OpenCL::Program* program, QObject* parent = nullptr);
   ::OpenCL::Event execute(
      ::OpenCL::CommandQueue* queue,
      int kernelSize,
      ::OpenCL::Buffer<Pairwise::Vector2>* in_data,
      cl_char clusterSize,
      ::OpenCL::Buffer<cl_char>* in_labels,
      cl_int sampleSize,
      cl_int minSamples,
      ::OpenCL::Buffer<cl_float>* out_correlations
   );
};



#endif
