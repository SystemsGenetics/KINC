#ifndef SIMILARITY_OPENCL_GMM_H
#define SIMILARITY_OPENCL_GMM_H
#include "similarity_opencl.h"



typedef struct
{
   cl_float pi;
   cl_float2 mu;
   cl_float4 sigma;
   cl_float4 sigmaInv;
   cl_float normalizer;
} cl_component;






/*!
 * This class implements the GMM kernel for the similarity analytic. This
 * kernel takes a list of pairwise data arrays and computes the number of
 * clusters and a list of cluster labels for each pair.
 */
class Similarity::OpenCL::GMM : public ::OpenCL::Kernel
{
   Q_OBJECT
public:
   /*!
    * Defines the arguments passed to the OpenCL kernel.
    */
   enum Argument
   {
      NumPairs
      ,Expressions
      ,SampleSize
      ,InIndex
      ,MinSamples
      ,MinClusters
      ,MaxClusters
      ,Criterion
      ,WorkX
      ,WorkN
      ,WorkLabels
      ,WorkComponents
      ,WorkMP
      ,WorkCounts
      ,WorkLogPi
      ,WorkGamma
      ,OutK
      ,OutLabels
   };
   explicit GMM(::OpenCL::Program* program, QObject* parent = nullptr);
   ::OpenCL::Event execute(
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
   );
};



#endif
