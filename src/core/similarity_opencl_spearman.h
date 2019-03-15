#ifndef SIMILARITY_OPENCL_SPEARMAN_H
#define SIMILARITY_OPENCL_SPEARMAN_H
#include "similarity_opencl.h"



/*!
 * This class implements the Pearson kernel for the similarity analytic. This
 * kernel takes a list of pairwise data arrays (with cluster labels) and computes
 * the Spearman correlation for each cluster in each pair.
 */
class Similarity::OpenCL::Spearman : public ::OpenCL::Kernel
{
   Q_OBJECT
public:
   /*!
    * Defines the arguments passed to the OpenCL kernel.
    */
   enum Argument
   {
      GlobalWorkSize
      ,InData
      ,ClusterSize
      ,InLabels
      ,SampleSize
      ,MinSamples
      ,WorkXY
      ,WorkRank
      ,OutCorrelations
   };
   explicit Spearman(::OpenCL::Program* program, QObject* parent = nullptr);
   ::OpenCL::Event execute(
      ::OpenCL::CommandQueue* queue,
      int globalWorkSize,
      int localWorkSize,
      ::OpenCL::Buffer<cl_float2>* in_data,
      cl_char clusterSize,
      ::OpenCL::Buffer<cl_char>* in_labels,
      cl_int sampleSize,
      cl_int minSamples,
      ::OpenCL::Buffer<cl_float>* work_xy,
      ::OpenCL::Buffer<cl_int>* work_rank,
      ::OpenCL::Buffer<cl_float>* out_correlations
   );
};



#endif
