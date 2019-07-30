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
      NumPairs
      ,Expressions
      ,SampleSize
      ,InIndex
      ,ClusterSize
      ,InLabels
      ,MinSamples
      ,WorkX
      ,WorkY
      ,OutCorrelations
   };
   explicit Spearman(::OpenCL::Program* program, QObject* parent = nullptr);
   ::OpenCL::Event execute(
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
   );
};



#endif
